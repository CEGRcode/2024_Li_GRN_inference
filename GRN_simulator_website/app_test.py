from flask import Flask, render_template, Response, request, jsonify, send_from_directory
import subprocess
import sys
import os
import zipfile
import shutil
import time, select
import logging
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import random
import math
import copy
import h5py
from scipy.integrate import solve_ivp
from distance_functions import *
from mutation_functions import *
from dynamics import *
from utility_functions import *
from population_functions import *
from GRN_Expanded_Combinatorial import GRN
from Combine_Redundant_Attractors import *
import io
import base64
from mpl_toolkits.mplot3d import Axes3D
import shlex
import json
import ast
from typing import List, Dict

matplotlib.use('Agg')

app = Flask(__name__)
input_value = None  # Store the input value globally

def unfold_edges_from_json(json_str: str) -> List[Dict]:
    """
    Parse a JSON string produced by your JS `JSON.stringify(...)` and
    return a list of Python dicts representing edges.
    Falls back to ast.literal_eval if json.loads fails (handles single quotes).
    """
    if json_str is None:
        return []
    json_str = json_str.strip()
    if json_str == "":
        return []

    try:
        data = json.loads(json_str)
        # JSON parsed successfully
    except json.JSONDecodeError:
        # Try a safer fallback for JS-like dict strings that use single quotes
        try:
            data = ast.literal_eval(json_str)
        except Exception as e:
            raise ValueError("Input is not valid JSON or Python literal") from e

    # Ensure we have a list
    if isinstance(data, dict):
        # If someone wrapped a single object, convert to list
        data = [data]
    if not isinstance(data, list):
        raise ValueError("Expected a list of edge objects")

    return data

def Json_grn_2_eq_string(grn_json):
    grn_data = unfold_edges_from_json(grn_json)
    regulation_content = {'independent': {'activation': [], 'inhibition': []}, 'synergistic': {'activation': {}, 'inhibition': {}}}
    for each in grn_data:
        if each['classes'][0] == 'solid':
            if each['classes'][1] == 'triangle':
                regulation_content['independent']['activation'].append(each['sourceData'])
            else:
                regulation_content['independent']['inhibition'].append(each['sourceData'])
        else:
            if each['classes'][1] == 'triangle':
                if each['label'] != '' and each['label'] not in regulation_content['synergistic']['activation']:
                    regulation_content['synergistic']['activation'][each['label']] = [each['sourceData']]
                elif each['label'] != '' and each['label'] in regulation_content['synergistic']['activation']:
                    regulation_content['synergistic']['activation'][each['label']].append(each['sourceData'])
                else:
                    pass
            else:
                if each['label'] != '' and each['label'] not in regulation_content['synergistic']['inhibition']:
                    regulation_content['synergistic']['inhibition'][each['label']] = [each['sourceData']]
                elif each['label'] != '' and each['label'] in regulation_content['synergistic']['inhibition']:
                    regulation_content['synergistic']['inhibition'][each['label']].append(each['sourceData'])
                else:
                    pass
    TempString_CSA_List = []
    TempString_CA1_sci = ''
    for each_act in regulation_content['independent']['activation']:
        TempString_CA1_sci = TempString_CA1_sci + 'modified_sigmoid(X, {}, {}, k, c1)*'.format(each_act['t1'], each_act['t2'])
        TempString_CA1_sci = TempString_CA1_sci[:-1]
        TempString_CSA_List.append(TempString_CA1_sci)
        TempString_CA1_sci = ''
    for each_act_complex in regulation_content['synergistic']['activation']:
        for each_act_subunit in regulation_content['synergistic']['activation'][each_act_complex]:
            TempString_CA1_sci = TempString_CA1_sci + 'modified_sigmoid(X, {}, {}, k, c1)*'.format(each_act_subunit['t1'], each_act_subunit['t2'])
        TempString_CA1_sci = TempString_CA1_sci[:-1]
        TempString_CSA_List.append(TempString_CA1_sci)
        TempString_CA1_sci = ''
    for each in TempString_CSA_List:
        TempString_CA1_sci = TempString_CA1_sci + '(1-{})*'.format(each)
    TempString_CA1_sci = TempString_CA1_sci[:-1]
    if TempString_CA1_sci == '':
        TempString_CA1_sci = '0'
    else:
        TempString_CA1_sci = '(1-' + TempString_CA1_sci + ')'
    TempString_CA2_sci = TempString_CA1_sci.replace('c1', 'c2')

    TempString_CSR_List = []
    TempString_CR1_sci = ''
    for each_inh in regulation_content['independent']['inhibition']:
        TempString_CR1_sci = TempString_CR1_sci + 'modified_sigmoid(Y, {}, {}, k, c3)*'.format(each_inh['t1'], each_inh['t2'])
        TempString_CR1_sci = TempString_CR1_sci[:-1]
        TempString_CSR_List.append(TempString_CR1_sci)
        TempString_CR1_sci = ''
    for each_inh_complex in regulation_content['synergistic']['inhibition']:
        for each_inh_subunit in regulation_content['synergistic']['inhibition'][each_inh_complex]:
            TempString_CR1_sci = TempString_CR1_sci + 'modified_sigmoid(Y, {}, {}, k, c3)*'.format(each_inh_subunit['t1'], each_inh_subunit['t2'])
        TempString_CR1_sci = TempString_CR1_sci[:-1]
        TempString_CSR_List.append(TempString_CR1_sci)
        TempString_CR1_sci = ''
    for each in TempString_CSR_List:
        TempString_CR1_sci = TempString_CR1_sci + '(1-{})*'.format(each)
    TempString_CR1_sci = TempString_CR1_sci[:-1]
    if TempString_CR1_sci == '':
        TempString_CR1_sci = '0'
    else:
        TempString_CR1_sci = '(1-' + TempString_CR1_sci + ')'
    TempString_CR2_sci = TempString_CR1_sci.replace('c3', 'c4')
    TempString_sci = 'f0+(1-f0)*{}-{}*f0+(f0p-1+f0)*{}*{}'.format(TempString_CA1_sci, TempString_CR1_sci, TempString_CA2_sci, TempString_CR2_sci)
    return TempString_sci

def f(X, Y, t1, t2, t3, k, f0, f0p, c1, c2, c3, c4):
    return (f0 +
            (1 - f0) * modified_sigmoid(X, t1, t2, k, c1) -
            modified_sigmoid(Y, t1, t3, k, c3) * f0 +
            (f0p - 1 + f0) * modified_sigmoid(X, t1, t2, k, c2) * modified_sigmoid(Y, t1, t3, k, c4))

@app.route('/')
def index():
    return render_template('GRN_Simulator_test.html')

@app.route('/unzip', methods=['POST'])
def unzip():
    data = request.get_json()
    filename = data.get('filename')
    if not filename:
        app.logger.error('No filename provided')
        return jsonify({'success': False, 'error': 'No filename provided'})

    file_path = os.path.join('static', 'data', filename)
    extract_path = os.path.join('static', 'data', 'unzipped')

    if not os.path.exists(file_path):
        app.logger.error('File not found: %s', file_path)
        return jsonify({'success': False, 'error': 'File not found'})

    try:
        with zipfile.ZipFile(file_path, 'r') as zip_ref:
            zip_ref.extractall(extract_path)
        return jsonify({'success': True})
    except Exception as e:
        app.logger.error('Error unzipping file: %s', e)
        return jsonify({'success': False, 'error': str(e)})

@app.route('/delete', methods=['POST'])
def delete_files():
    time.sleep(15)
    extract_path = os.path.join('static', 'data', 'unzipped')
    try:
        if os.path.exists(extract_path):
            shutil.rmtree(extract_path)
        return jsonify({'success': True})
    except Exception as e:
        app.logger.error('Error deleting files: %s', e)
        return jsonify({'success': False, 'error': str(e)})
    

@app.route('/set-input', methods=['POST'])
def set_input():
    global input_value
    input_value = request.form['input_value']
    return '', 204  # No Content

@app.route('/api/plot', methods=['POST'])
def plot_api():
    try:
        params = request.json
        html = generate_plot(**params)
        return jsonify({
            'success': True,
            'html': html,
            'params': params
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

@app.route('/stream-data')
def stream_data():
    global input_value

    def generate():
        # Launch the subprocess in unbuffered mode (-u) so that output is sent immediately.
        with subprocess.Popen(['python', '-u', 'GRN_simulator.py', input_value],
                              stdout=subprocess.PIPE, text=True) as process:
            last_heartbeat = time.time()
            
            # Read output in a loop, non-blockingly
            while True:
                # Wait for up to 0.1 seconds for data on process.stdout
                reads, _, _ = select.select([process.stdout], [], [], 0.1)
                
                if reads:
                    line = process.stdout.readline()
                    if not line:  # EOF reached
                        break
                    if line.strip() == 'Process completed.':
                        yield f"data: {line}\n\n"  # Signal completion
                        break
                    yield f"data: {line}\n\n"
                    last_heartbeat = time.time()  # Reset heartbeat timer on new output
                else:
                    # If no new output and more than 2 seconds have passed, send a heartbeat.
                    if time.time() - last_heartbeat > 2:
                        yield "data: heartbeat\n\n"
                        last_heartbeat = time.time()
                
                # Break out of the loop if the process has terminated.
                if process.poll() is not None:
                    break

            # Process any remaining output lines.
            for line in process.stdout:
                yield f"data: {line}\n\n"
            process.stdout.close()
            process.wait()
        yield "data: Attractor calculation completed.\n\n"

    return Response(generate(), mimetype='text/event-stream')

@app.route('/runtime', methods=['POST'])
def runtime():
    input_data = request.get_json()['input']
    
    # Run the command
    try:
        result = subprocess.run(['python', 'GRN_simulator.py', input_data], capture_output=True, text=True, check=False)
        output = result.stdout
        error = result.stderr
    except subprocess.CalledProcessError as e:
        output = e.stdout
        error = e.stderr
    
    return jsonify({
        'output': output,
        'error': error
    })

@app.route('/composite_input', methods=['POST'])
def composite_input():
    # Retrieve the data from the request
    data = request.form.get('input_samples')
    if data is None:
        return jsonify({'error': 'No input_samples provided'}), 400
    try:
        # Split and process input_samples
        source, target = data.split('\t')
        Selected_Sample_Codes = [source]
        target = target.split(',')

        # Convert the input values to a list of floats
        for each in target:
            Selected_Sample_Codes.append(each)
        # Create the plot
        # Generate 10 distinct colors
        plt.figure(figsize=(40, 20))
        colors = generate_distinct_colors(10)
        YEP_GO = {}
        YEP_GO_file = open('./static/data/YEP_GO.txt', 'r')
        for line in YEP_GO_file:
            YEP_GO[line.split('\t')[0].upper()] = line.split('\t')[1][:-1]
        YEP_GO_file.close()

        SampleKeyFile = open('./static/data/SupplementaryData-Table4_Sample-Key_tabular.tab', 'r')
        Sample_Keys = {}
        for line in SampleKeyFile:
            Sample_Keys[int(line.split()[0])] = line.split()[1]
        SampleKeyFile.close()

        # obtain the binding sites intersect
        Binding_sites = {}
        infile = open('./static/data/Binding_sites.txt', 'r')
        for line in infile:
            Binding_sites[line.split()[0]] = line.split()[1]
        infile.close()

        composite_command = 'python plotter.py'        
        BED_ID = Selected_Sample_Codes[0]
        colors = generate_distinct_colors(len(Selected_Sample_Codes))
        forward_bg, reverse_bg = read_cdt_h5(BED_ID, BED_ID)
        with open('./static/data/temp_composite_ref.out', 'w') as outfile:
            outfile.write('\t'+'\t'.join(map(str, [tick_ for tick_ in range(-250, 251)]))+'\n')
            outfile.write(BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)])+'_sense'+'\t'+'\t'.join(map(str, forward_bg))+'\n')
            outfile.write(BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)])+'_anti'+'\t'+'\t'.join(map(str, reverse_bg))+'\n')
        composite_command = composite_command + ' composite ./static/data/temp_composite_ref.out --color \'#C4A77A\' --name {}'.format(BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)]))
        i = 0
        for each in Selected_Sample_Codes[1:]:
            forward_CDT, reverse_CDT = read_cdt_h5(BED_ID, each)
            color_i = "#{:02X}{:02X}{:02X}".format(int(round(colors[i][0]*255)), int(round(colors[i][1]*255)), int(round(colors[i][2]*255)))
            with open('./static/data/temp_composite_{}.out'.format(i), 'w') as outfile:
                outfile.write('\t'+'\t'.join(map(str, [tick_ for tick_ in range(-250, 251)]))+'\n')
                outfile.write(each+'_{}'.format(Sample_Keys[int(each)])+'_sense'+'\t'+'\t'.join(map(str, forward_CDT))+'\n')
                outfile.write(each+'_{}'.format(Sample_Keys[int(each)])+'_anti'+'\t'+'\t'.join(map(str, reverse_CDT))+'\n')
            composite_command = composite_command + ' composite ./static/data/temp_composite_{}.out --color \'{}\' --name {}'.format(i, color_i, each+'_{}'.format(Sample_Keys[int(each)]))
            i=i+1
        forward_IgG, reverse_IgG = read_cdt_h5(BED_ID, 'IgG')
        IgG_scaling_factor = 9 * max(max(forward_IgG), abs(min(reverse_IgG))) / max(max(forward_bg), abs(min(reverse_bg)))
        forward_IgG = list(np.array(forward_IgG) / IgG_scaling_factor)
        reverse_IgG = list(np.array(reverse_IgG) / IgG_scaling_factor)
        with open('./static/data/temp_composite_igg.out', 'w') as outfile:
            outfile.write('\t'+'\t'.join(map(str, [tick_ for tick_ in range(-250, 251)]))+'\n')
            outfile.write(BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)])+'_sense'+'\t'+'\t'.join(map(str, forward_IgG))+'\n')
            outfile.write(BED_ID+'_{}'.format(Sample_Keys[int(BED_ID)])+'_anti'+'\t'+'\t'.join(map(str, reverse_IgG))+'\n')
        composite_command = composite_command + ' composite ./static/data/temp_composite_igg.out --color \'#bfbfbb\' --name {}'.format('IgG')
        if Sample_Keys[int(BED_ID)].capitalize()+'_'+Sample_Keys[int(each)].capitalize() not in Binding_sites:
            plot_xlabel = '{} bound sites for {}, of which {} are bound by {}'.format(Binding_sites[Sample_Keys[int(BED_ID)].capitalize()], Sample_Keys[int(BED_ID)].capitalize(), 'none', Sample_Keys[int(each)].capitalize())
        else:
            plot_xlabel = '{} bound sites for {}, of which {} are bound by {}'.format(Binding_sites[Sample_Keys[int(BED_ID)].capitalize()], Sample_Keys[int(BED_ID)].capitalize(), Binding_sites[Sample_Keys[int(BED_ID)].capitalize()+'_'+Sample_Keys[int(each)].capitalize()], Sample_Keys[int(each)].capitalize())
        composite_command = composite_command + ' plot --xlabel \"' + plot_xlabel + '\"' + ' --title \"' + 'Overlapping composite plot on {} bound sites'.format(Sample_Keys[int(BED_ID)].capitalize()) + '\" --smoothing 3 --color-trace --out ./static/data/dashed_{}.svg'.format(str(source)+'_'+str(target[0]))        
        composite_args = shlex.split(composite_command)
        proc = subprocess.run(composite_args, check=False, capture_output=True, text=True)
        # Example processing
        result = f"{source}_{target[0]}"

        # Return a JSON response with the result
        return jsonify({'message': result}), 200
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
# Route to read the ChIP-exo sample index file content
BASE_DIR = os.path.abspath(os.path.dirname(__file__))
FILE_PATH = BASE_DIR + '/static/data/YEP_best_rep.txt'
@app.route('/read-file', methods=['GET'])
def read_file():
    try:
        with open(FILE_PATH, 'r') as file:
            content = file.read()
        return jsonify({'content': content})
    except FileNotFoundError:
        return jsonify({'error': 'File not found'}), 404

@app.route('/plot', methods=['GET'])
def plot():
    try:
        node_name = str(request.args.get('label', 'factor X'))
        t1   = float(request.args.get('t1', 0.33))
        t2   = float(request.args.get('t2', 0.66))
        t3   = float(request.args.get('t3', 0.66))
        k    = float(request.args.get('k', 35))
        f0   = float(request.args.get('f0', 0.5))
        f0p  = float(request.args.get('f0p', 0.5))
        c1   = float(request.args.get('c1', 1))
        c2   = float(request.args.get('c2', 1))
        c3   = float(request.args.get('c3', 1))
        c4   = float(request.args.get('c4', 1))
        elev = float(request.args.get('elev', 25))
        azim = float(request.args.get('azim', 160))
        grn = request.args.get('grn', '')
    except Exception as e:
        return jsonify({'error': str(e)})

    X, Y = np.meshgrid(np.linspace(0, 1, 100), np.linspace(0, 1, 100))
    y = eval(Json_grn_2_eq_string(grn))
    if type(y) == float:
        y = np.full_like(X, y)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, y, cmap='viridis', edgecolor='none')
    ax.view_init(elev=elev, azim=azim)

    ax.set_xlabel('Comb. Hill output of activator(s)')
    ax.set_ylabel('Comb. Hill output of inhibitor(s)')
    ax.set_zlabel('TG expression')

    ax.set_title('Regulation function of {}'.format(node_name))
    
    img = io.BytesIO()
    plt.savefig(img, format='png', dpi=300)
    plt.close(fig)
    img.seek(0)
    
    img_b64 = base64.b64encode(img.getvalue()).decode('utf8')
    return jsonify({'img': img_b64})

if __name__ == '__main__':
    app.run(
        ssl_context=('/etc/letsencrypt/live/grn.cac.cornell.edu/fullchain.pem',
        '/etc/letsencrypt/live/grn.cac.cornell.edu/privkey.pem'),
        host='0.0.0.0',
        port=5000
    )
