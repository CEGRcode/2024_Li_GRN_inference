from flask import Flask, render_template, request, jsonify
import subprocess
import sys
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import random
import math
import copy
from scipy.integrate import solve_ivp
from distance_functions import *
from mutation_functions import *
from dynamics import *
from utility_functions import *
from population_functions import *
from GRN_Expanded_Combinatorial import GRN
from Combine_Redundant_Attractors import *

app = Flask(__name__, static_folder='/home/ubuntu/GRN_UI_app_test/public')

@app.route('/')
def index():
    return render_template('GRN_Simulator_test.html')

@app.route('/run_python', methods=['POST'])
def run_python():
    input_data = request.json.get('input_data')
    
    # Path to the Python executable
    python_script = 'GRN_simulator_test.py'  # Change this to the path of your script
    
    try:
        # Pass input_data as an argument to the Python script
        result = subprocess.run(['python', python_script, input_data], capture_output=True, text=True, check=True)
        output = result.stdout
    except subprocess.CalledProcessError as e:
        output = f"An error occurred: {e.output}"
    
    return jsonify({'output': output})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
