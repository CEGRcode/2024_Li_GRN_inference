import copy
from scipy.integrate import solve_ivp
import numpy as np
import sympy as sp
from scipy.optimize import minimize
import random

def ActivatorSigmoid(x, t, k):
    '''Activating Hill function: threshold t, Hill coefficient k'''
    return (x**k / (x**k + t**k))

def modified_sigmoid(x, t1_, t2_, k, c):
    # Two Hill functions with different thresholds t1 and t2
    sigma1 = ActivatorSigmoid(x, t1_, k)
    sigma2 = ActivatorSigmoid(x, t2_, k)
    return (sigma1 + c*sigma2) / (1+c)

def GetAttractorDistance(mRNA1, WT_mRNA, transcriptionalprofilemax):
    '''Compute normalized distance between an two attractors'''
    outdistance = 0.0
    for y in range(0, len(WT_mRNA)):
        outdistance = outdistance + \
            abs((mRNA1[y]-WT_mRNA[y])/transcriptionalprofilemax[y])
    return outdistance

def Run_Dynamics(GRN_instance_ori, index_i_func, index_j_func, TrainningCount_func, WTTP, TRMAX, return_dict):
    '''Simulate dynamics'''
    t_ticks = [int(TrainningCount_func/250)*tick for tick in range(0, 250+1)]
    GRN_instance = copy.deepcopy(GRN_instance_ori)
    specific_gene = int(GRN_instance.Name.split('_')[1])
    GRN_instance.SetmRNA(WTTP[str(index_i_func)][1], WTTP[str(index_i_func)][0], TRMAX)
    scipystring = GRN_instance.Delta_mRNA(WTTP[str(index_i_func)][0], TRMAX, list(GRN_instance.mRNA))
    scipystring = (scipystring
                   + '\nInitialState = [GRN_instance.mRNA[{}]]'.format(specific_gene)
                   + '\nsol = solve_ivp(update_mRNA_protein, [0, TrainningCount_func], InitialState, '
                   +                    'args=([WTTP[str(index_i_func)][0]]), t_eval=t_ticks, method=\'RK45\')')
    scipystring = scipystring + '\nGRN_instance.sol = sol.y'
    #print(scipystring, '\n')
    exec(scipystring)
    sol = GRN_instance.sol
    GRN_instance.sol = []
    sol = [int(x.item()) for x in list(np.round(sol.T, 2))]
    #print('sol: ', sol)
    return_dict[index_j_func] = sol
    sol = []
    GRN_instance = []
    return

def get_param_by_TPM(TPMs):
    '''Obtain the f0, f0p, c1, c2, c3, c4, and the threshold ts for each regulated gene.'''

    max_number = sorted(TPMs)[-1]
    min_number = sorted(TPMs)[0]
    rest_numbers = sorted(TPMs)[1:-1] / max_number # excluded the maximal and minimal which are fixed as 1 and 0 in the function.
    model_expression_levels = ['f0/(1+c3)', '((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))', 'f0p', 'f0',
                               'f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))',
                               '(c1+f0)/(c1+1)', '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)']
    initial_guess = [0.5, 0.5, 1, 1, 1, 1]
    bounds = [(0, 1), (0, 1), (0, None), (0, None), (0, None), (0, None)]

    if len(rest_numbers) == 6:
        model_expression_levels = ['f0/(1+c3)', '((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))', 'f0p', 'f0',
                                   'f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))',
                                   '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)', '(c1+f0)/(c1+1)']
        to_remove = [['f0/(1+c3)'], ['((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))'], ['(c1+f0)/(c1+1)'],
                     ['1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)']]
        scores_ = []
        constrain_string_ = [
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0p-f0/(1+c3)',
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0-((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))',
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0p-(c1+f0)/(c1+1)',
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0-1+(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)'
        ]
        eq0_string_ = ['eq0 = sp.Eq(f0p-f0/(1+c3), 0)',
                       'eq0 = sp.Eq(f0-((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1))), 0)',
                       'eq0 = sp.Eq(f0p-(c1+f0)/(c1+1), 0)',
                       'eq0 = sp.Eq(f0-1+(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1), 0)']
        results_list = []
        for to_remove_i in range(0, len(to_remove)):
            constrain_string = constrain_string_[to_remove_i]
            eq0_string = eq0_string_[to_remove_i]
            model_expression_levels_ = [item for item in model_expression_levels if item not in to_remove[to_remove_i]]
            least_square_string = 'def objective(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n'
            return_string = '\treturn '
            f0, f0p, c1, c2, c3, c4 = sp.symbols('f0, f0p, c1, c2, c3, c4')
            for i in range(0, len(rest_numbers)):
                least_square_string = least_square_string + '\tr{} = {} - {}\n'.format(i + 1,
                                                                                       model_expression_levels_[i],
                                                                                       rest_numbers[i])
                return_string = return_string + 'r{}**2+'.format(i + 1)
                exec('eq{} = sp.Eq({}, {})'.format(i + 1, model_expression_levels_[i], rest_numbers[i]))
            exec(eq0_string)
            return_string = return_string[:-1]
            least_square_string = least_square_string + return_string
            # print(eq00, '\n', eq0, '\n', eq1, '\n', eq2, '\n', eq3, '\n', eq4, '\n', eq5)
            try:
                solutions = sp.solve((eq00, eq0, eq1, eq2, eq3, eq4, eq5), (f0, f0p, c1, c2, c3, c4))
                print('ori solution ==>: ', solutions)
                solutions = list(solutions[0])
                values = {f0: 0.5, f0p: 0.5, c1: 1, c2: 1, c3: 1, c4: 1}
                result_f0_f0p_c1_c2_c3_c4 = [
                    term if isinstance(term, (sp.core.numbers.Float)) else float(term.subs(values)) for term in
                    solutions]
            except:
                solutions = None
            if solutions == None:
                exec(least_square_string, globals())
                exec(constrain_string, globals())
                cons = [{'type': 'eq', 'fun': constraint}]
                result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
                results_list.append(result.x)
                examination = [result.x[0] / (1 + result.x[4]),
                               ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                           (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                               result.x[1],
                               result.x[0],
                               result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                           result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                           result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                           result.x[5] / (1 + result.x[5])),
                               (result.x[2] + result.x[0]) / (result.x[2] + 1),
                               1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (
                                           result.x[5] / (result.x[5] + 1)) * (result.x[1] + result.x[0] - 1)]
                if min(examination) < 0 or max(examination) > 1:
                    scores_.append(1e6)
                else:
                    scores_.append(result.fun)
            elif min(result_f0_f0p_c1_c2_c3_c4) >= 0 and max(result_f0_f0p_c1_c2_c3_c4[0],
                                                             result_f0_f0p_c1_c2_c3_c4[1]) <= 1:
                results_list.append(result_f0_f0p_c1_c2_c3_c4)
                scores_.append(0)
            else:
                exec(least_square_string, globals())
                exec(constrain_string, globals())
                cons = [{'type': 'eq', 'fun': constraint}]
                result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
                results_list.append(result.x)
                examination = [result.x[0] / (1 + result.x[4]),
                               ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                           (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                               result.x[1],
                               result.x[0],
                               result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                           result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                           result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                           result.x[5] / (1 + result.x[5])),
                               (result.x[2] + result.x[0]) / (result.x[2] + 1),
                               1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (
                                           result.x[5] / (result.x[5] + 1)) * (result.x[1] + result.x[0] - 1)]
                if min(examination) < 0 or max(examination) > 1:
                    scores_.append(1e6)
                else:
                    scores_.append(result.fun)
            # print(least_square_string)
        min_score_index = np.argmin(scores_)
        result_t = [0.66, 0.33, 0.33]
        result_f0_f0p_c1_c2_c3_c4 = results_list[min_score_index]

    elif len(rest_numbers) == 5:
        model_expression_levels = ['f0/(1+c3)', '((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))', 'f0p', 'f0',
                                   'f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))',
                                   '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)', '(c1+f0)/(c1+1)']
        to_remove = [['f0', 'f0/(1+c3)'], ['f0', '(c1+f0)/(c1+1)'],
                     ['f0p', '((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))'],
                     ['f0p', '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)']]
        scores_ = []
        constrain_string_ = [
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0',
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn 1-f0',
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0p',
            'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn 1-f0p'
        ]
        constrain_string_additional_ = [
            'def constraint_2(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn f0/(1+c3)',
            'def constraint_2(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn 1-(c1+f0)/(c1+1)',
            'def constraint_2(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn ((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))',
            'def constraint_2(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn (f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)'
        ]
        eq0_string_ = ['eq0 = sp.Eq(f0, 0)', 'eq0 = sp.Eq(f0, 1)', 'eq0 = sp.Eq(f0p, 0)', 'eq0 = sp.Eq(f0p, 1)']
        eq00_string_ = ['eq00 = sp.Eq(f0/(1+c3), 0)', 'eq00 = sp.Eq((c1+f0)/(c1+1), 1)',
                        'eq00 = sp.Eq(((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1))), 0)',
                        'eq00 = sp.Eq(1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1), 1)']
        results_list = []
        for to_remove_i in range(0, len(to_remove)):
            constrain_string = constrain_string_[to_remove_i]
            constrain_string_additional = constrain_string_additional_[to_remove_i]
            eq0_string = eq0_string_[to_remove_i]
            eq00_string = eq00_string_[to_remove_i]
            model_expression_levels_ = [item for item in model_expression_levels if item not in to_remove[to_remove_i]]
            least_square_string = 'def objective(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n'
            return_string = '\treturn '
            f0, f0p, c1, c2, c3, c4 = sp.symbols('f0, f0p, c1, c2, c3, c4')
            for i in range(0, len(rest_numbers)):
                least_square_string = least_square_string + '\tr{} = {} - {}\n'.format(i + 1,
                                                                                       model_expression_levels_[i],
                                                                                       rest_numbers[i])
                return_string = return_string + 'r{}**2+'.format(i + 1)
                exec('eq{} = sp.Eq({}, {})'.format(i + 1, model_expression_levels_[i], rest_numbers[i]))
            exec(eq0_string)
            exec(eq00_string)
            return_string = return_string[:-1]
            least_square_string = least_square_string + return_string
            # print(eq00, '\n', eq0, '\n', eq1, '\n', eq2, '\n', eq3, '\n', eq4, '\n', eq5)
            try:
                solutions = sp.solve((eq00, eq0, eq1, eq2, eq3, eq4, eq5), (f0, f0p, c1, c2, c3, c4))
                print('ori solution ==>: ', solutions)
                solutions = list(solutions[0])
                values = {f0: 0.5, f0p: 0.5, c1: 1, c2: 1, c3: 1, c4: 1}
                result_f0_f0p_c1_c2_c3_c4 = [
                    term if isinstance(term, (sp.core.numbers.Float)) else float(term.subs(values)) for term in
                    solutions]
            except:
                solutions = None
            if solutions == None:
                exec(least_square_string, globals())
                exec(constrain_string, globals())
                exec(constrain_string_additional, globals())
                cons = [{'type': 'eq', 'fun': constraint}, {'type': 'eq', 'fun': constraint_2}]
                result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
                results_list.append(result.x)
                examination = [result.x[0] / (1 + result.x[4]),
                               ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                           (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                               result.x[1],
                               result.x[0],
                               result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                           result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                           result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                           result.x[5] / (1 + result.x[5])),
                               (result.x[2] + result.x[0]) / (result.x[2] + 1),
                               1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (
                                           result.x[5] / (result.x[5] + 1)) * (result.x[1] + result.x[0] - 1)]
                if min(examination) < 0 or max(examination) > 1:
                    scores_.append(1e6)
                else:
                    scores_.append(result.fun)
            elif min(result_f0_f0p_c1_c2_c3_c4) >= 0 and max(result_f0_f0p_c1_c2_c3_c4[0],
                                                             result_f0_f0p_c1_c2_c3_c4[1]) <= 1:
                results_list.append(result_f0_f0p_c1_c2_c3_c4)
                scores_.append(0)
            else:
                exec(least_square_string, globals())
                exec(constrain_string, globals())
                exec(constrain_string_additional, globals())
                cons = [{'type': 'eq', 'fun': constraint}, {'type': 'eq', 'fun': constraint_2}]
                result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
                results_list.append(result.x)
                examination = [result.x[0] / (1 + result.x[4]),
                               ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                           (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                               result.x[1],
                               result.x[0],
                               result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                           result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                           result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                           result.x[5] / (1 + result.x[5])),
                               (result.x[2] + result.x[0]) / (result.x[2] + 1),
                               1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (
                                           result.x[5] / (result.x[5] + 1)) * (result.x[1] + result.x[0] - 1)]
                if min(examination) < 0 or max(examination) > 1:
                    scores_.append(1e6)
                else:
                    scores_.append(result.fun)
            # print(least_square_string)

        min_score_index = np.argmin(scores_)
        result_t = [0.66, 0.33, 0.33]
        result_f0_f0p_c1_c2_c3_c4 = results_list[min_score_index]

    elif len(rest_numbers) == 4:
        result_t = [0.66, 0.33, 0.66]
        to_remove = ['f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))', 'f0/(1+c3)',
                     '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)']
        model_expression_levels = ['((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))', 'f0p', 'f0', '(c1 + f0)/(c1 + 1)']
        least_square_string = 'def objective(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n'
        return_string = '\treturn '
        f0, f0p, c1, c2, c3, c4 = sp.symbols('f0, f0p, c1, c2, c3, c4')
        for i in range(0, len(rest_numbers)):
            least_square_string = least_square_string + '\tr{} = {} - {}\n'.format(i + 1, model_expression_levels[i],
                                                                                   rest_numbers[i])
            return_string = return_string + 'r{}**2+'.format(i + 1)
            exec('eq{} = sp.Eq({}, {})'.format(i + 1, model_expression_levels[i], rest_numbers[i]))
        return_string = return_string[:-1]
        least_square_string = least_square_string + return_string
        exec(least_square_string, globals())
        result_f0_f0p_c1_c2_c3_c4 = []
        try:
            solutions = sp.solve((eq1, eq2, eq3, eq4), (f0, f0p, c1, c2, c3, c4))
            solutions = list(solutions[0])
            values = {f0: 0.5, f0p: 0.5, c1: 1, c2: 1, c3: 1, c4: 1}
            result_f0_f0p_c1_c2_c3_c4 = [term if isinstance(term, (sp.core.numbers.Float)) else float(term.subs(values))
                                         for term in solutions]
        except:
            solutions = None
        if solutions == None:
            exec(least_square_string, globals())
            result = minimize(objective, initial_guess, bounds=bounds, method='SLSQP')
            result_list_1 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_1 = 1e6
            else:
                scores_1 = result.fun
        elif min(result_f0_f0p_c1_c2_c3_c4) >= 0 and max(result_f0_f0p_c1_c2_c3_c4[0],
                                                         result_f0_f0p_c1_c2_c3_c4[1]) <= 1:
            result_list_1 = result_f0_f0p_c1_c2_c3_c4
            scores_1 = 0
        else:
            exec(least_square_string, globals())
            result = minimize(objective, initial_guess, bounds=bounds, method='SLSQP')
            result_list_1 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_1 = 1e6
            else:
                scores_1 = result.fun
        result_t = [0.66, 0.66, 0.33]
        to_remove = ['f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))', '(c1+f0)/(c1+1)',
                     '((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))']
        model_expression_levels = ['f0/(1+c3)', 'f0', 'f0p', '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)']
        least_square_string = 'def objective(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n'
        return_string = '\treturn '
        f0, f0p, c1, c2, c3, c4 = sp.symbols('f0, f0p, c1, c2, c3, c4')
        for i in range(0, len(rest_numbers)):
            least_square_string = least_square_string + '\tr{} = {} - {}\n'.format(i + 1, model_expression_levels[i],
                                                                                   rest_numbers[i])
            return_string = return_string + 'r{}**2+'.format(i + 1)
            exec('eq{} = sp.Eq({}, {})'.format(i + 1, model_expression_levels[i], rest_numbers[i]))
        return_string = return_string[:-1]
        least_square_string = least_square_string + return_string
        exec(least_square_string, globals())
        result_f0_f0p_c1_c2_c3_c4 = []
        try:
            solutions = sp.solve((eq1, eq2, eq3, eq4), (f0, f0p, c1, c2, c3, c4))
            solutions = list(solutions[0])
            values = {f0: 0.5, f0p: 0.5, c1: 1, c2: 1, c3: 1, c4: 1}
            result_f0_f0p_c1_c2_c3_c4 = [term if isinstance(term, (sp.core.numbers.Float)) else float(term.subs(values))
                                         for term in solutions]
        except:
            solutions = None
        if solutions == None:
            exec(least_square_string, globals())
            result = minimize(objective, initial_guess, bounds=bounds, method='SLSQP')
            result_list_2 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_2 = 1e6
            else:
                scores_2 = result.fun
        elif min(result_f0_f0p_c1_c2_c3_c4) >= 0 and max(result_f0_f0p_c1_c2_c3_c4[0],
                                                         result_f0_f0p_c1_c2_c3_c4[1]) <= 1:
            result_list_2 = result_f0_f0p_c1_c2_c3_c4
            scores_2 = 0
        else:
            exec(least_square_string, globals())
            result = minimize(objective, initial_guess, bounds=bounds, method='SLSQP')
            result_list_2 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_2 = 1e6
            else:
                scores_2 = result.fun
        if scores_1 <= scores_2:
            result_f0_f0p_c1_c2_c3_c4 = result_list_1
            result_t = [0.66, 0.33, 0.66]
        else:
            result_f0_f0p_c1_c2_c3_c4 = result_list_2
            result_t = [0.66, 0.66, 0.33]

    elif len(rest_numbers) == 3:
        result_t = [0.66, 0.66, 0.33]
        to_remove = ['1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)',
                     'f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))', '(c1+f0)/(c1+1)',
                     '((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))']
        constrain_string = 'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn 1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)-f0'
        eq0_string = 'eq0 = sp.Eq(1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)-f0, 0)'
        model_expression_levels = ['f0p', 'f0/(1+c3)', 'f0']
        least_square_string = 'def objective(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n'
        return_string = '\treturn '
        f0, f0p, c1, c2, c3, c4 = sp.symbols('f0, f0p, c1, c2, c3, c4')
        for i in range(0, len(rest_numbers)):
            least_square_string = least_square_string + '\tr{} = {} - {}\n'.format(i + 1, model_expression_levels[i],
                                                                                   rest_numbers[i])
            return_string = return_string + 'r{}**2+'.format(i + 1)
            exec('eq{} = sp.Eq({}, {})'.format(i + 1, model_expression_levels[i], rest_numbers[i]))
        exec(eq0_string)
        exec(constrain_string, globals())
        return_string = return_string[:-1]
        least_square_string = least_square_string + return_string
        # print(eq0, '\n', eq1, '\n', eq2, '\n', eq3, '\n')
        result_f0_f0p_c1_c2_c3_c4 = []
        try:
            solutions = sp.solve((eq0, eq1, eq2, eq3), (f0, f0p, c1, c2, c3, c4))
            solutions = list(solutions[0])
            values = {f0: 0.5, f0p: 0.5, c1: 1, c2: 1, c3: 1, c4: 1}
            result_f0_f0p_c1_c2_c3_c4 = [term if isinstance(term, (sp.core.numbers.Float)) else float(term.subs(values))
                                         for term in solutions]
        except:
            solutions = None
        if solutions == None:
            exec(least_square_string, globals())
            cons = [{'type': 'eq', 'fun': constraint}]
            result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
            result_f0_f0p_c1_c2_c3_c4 = result.x
            result_list_1 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_1 = 1e6
            else:
                scores_1 = result.fun
        elif min(result_f0_f0p_c1_c2_c3_c4) >= 0 and max(result_f0_f0p_c1_c2_c3_c4[0],
                                                         result_f0_f0p_c1_c2_c3_c4[1]) <= 1:
            scores_1 = 0
            result_list_1 = result_f0_f0p_c1_c2_c3_c4
            pass
        else:
            exec(least_square_string, globals())
            cons = [{'type': 'eq', 'fun': constraint}]
            result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
            result_f0_f0p_c1_c2_c3_c4 = result.x
            result_list_1 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_1 = 1e6
            else:
                scores_1 = result.fun

        result_t = [0.66, 0.33, 0.66]
        to_remove = ['(c1+f0)/(c1+1)', 'f0+((1-f0)*(c1/(c1+1)))-(f0*c3/(1+c3))+(f0+f0p-1)*(c2/(1+c2))*(c4/(1+c4))',
                     'f0/(1+c3)', '1-(f0*c3/(c3+1))+(c4/(c4+1))*(f0p+f0-1)']
        constrain_string = 'def constraint(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n\treturn c1+f0 - f0p*(c1+1)'
        eq0_string = 'eq0 = sp.Eq((c1 + f0)/(c1 + 1)-f0p, 0)'
        model_expression_levels = ['((1-f0)*(c1/(1+c1)))+((f0p-1+f0)*(c2/(c2+1)))', 'f0', 'f0p']
        least_square_string = 'def objective(params):\n\tf0, f0p, c1, c2, c3, c4 = params\n'
        return_string = '\treturn '
        f0, f0p, c1, c2, c3, c4 = sp.symbols('f0, f0p, c1, c2, c3, c4')
        for i in range(0, len(rest_numbers)):
            least_square_string = least_square_string + '\tr{} = {} - {}\n'.format(i + 1, model_expression_levels[i],
                                                                                   rest_numbers[i])
            return_string = return_string + 'r{}**2+'.format(i + 1)
            exec('eq{} = sp.Eq({}, {})'.format(i + 1, model_expression_levels[i], rest_numbers[i]))
        exec(eq0_string)
        exec(constrain_string, globals())
        return_string = return_string[:-1]
        least_square_string = least_square_string + return_string
        # print(eq0, '\n', eq1, '\n', eq2, '\n', eq3, '\n')
        result_f0_f0p_c1_c2_c3_c4 = []
        try:
            solutions = sp.solve((eq0, eq1, eq2, eq3), (f0, f0p, c1, c2, c3, c4))
            solutions = list(solutions[0])
            values = {f0: 0.5, f0p: 0.5, c1: 1, c2: 1, c3: 1, c4: 1}
            result_f0_f0p_c1_c2_c3_c4 = [term if isinstance(term, (sp.core.numbers.Float)) else float(term.subs(values))
                                         for term in solutions]
        except:
            solutions = None
        if solutions == None:
            exec(least_square_string, globals())
            cons = [{'type': 'eq', 'fun': constraint}]
            result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
            result_f0_f0p_c1_c2_c3_c4 = result.x
            result_list_2 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_2 = 1e6
            else:
                scores_2 = result.fun
        elif min(result_f0_f0p_c1_c2_c3_c4) >= 0 and max(result_f0_f0p_c1_c2_c3_c4[0],
                                                         result_f0_f0p_c1_c2_c3_c4[1]) <= 1:
            scores_2 = 0
            result_list_2 = result_f0_f0p_c1_c2_c3_c4
            pass
        else:
            exec(least_square_string, globals())
            cons = [{'type': 'eq', 'fun': constraint}]
            result = minimize(objective, initial_guess, bounds=bounds, constraints=cons, method='SLSQP')
            result_f0_f0p_c1_c2_c3_c4 = result.x
            result_list_2 = result.x
            examination = [result.x[0] / (1 + result.x[4]),
                           ((1 - result.x[0]) * (result.x[2] / (1 + result.x[2]))) + (
                                       (result.x[1] - 1 + result.x[0]) * (result.x[3] / (result.x[3] + 1))),
                           result.x[1],
                           result.x[0],
                           result.x[0] + ((1 - result.x[0]) * (result.x[2] / (result.x[2] + 1))) - (
                                       result.x[0] * result.x[4] / (1 + result.x[4])) + (
                                       result.x[0] + result.x[1] - 1) * (result.x[3] / (1 + result.x[3])) * (
                                       result.x[5] / (1 + result.x[5])),
                           (result.x[2] + result.x[0]) / (result.x[2] + 1),
                           1 - (result.x[0] * result.x[4] / (result.x[4] + 1)) + (result.x[5] / (result.x[5] + 1)) * (
                                       result.x[1] + result.x[0] - 1)]
            if min(examination) < 0 or max(examination) > 1:
                scores_2 = 1e6
            else:
                scores_2 = result.fun
        if scores_1 <= scores_2:
            result_f0_f0p_c1_c2_c3_c4 = result_list_1
            result_t = [0.66, 0.66, 0.33]
        else:
            result_f0_f0p_c1_c2_c3_c4 = result_list_2
            result_t = [0.66, 0.33, 0.66]

    elif len(rest_numbers) == 2:
        model_expression_levels = ['f0', 'f0p']
        result_t = [0.5, 0.5, 0.5]
        result_f0_f0p_c1_c2_c3_c4 = [rest_numbers[0], rest_numbers[1], 1, 1, 1, 1]

    elif len(rest_numbers) == 1:
        result_t = [0.5, 0.5, 0.5]
        result_f0_f0p_c1_c2_c3_c4 = [rest_numbers[0], rest_numbers[0], 1, 1, 1, 1]

    elif len(rest_numbers) == 0:
        result_t = [0.5, 0.5, 0.5]
        result_f0_f0p_c1_c2_c3_c4 = random.choice([[0, 0, 1, 1, 1, 1], [1, 1, 1, 1, 1, 1]])
    else:
        pass

    return (result_t, list(result_f0_f0p_c1_c2_c3_c4))

def generate_plot(f0=0.5, t=0.5, k=25, elev=25, azim=160):
    """Make 3-D plot for interactive parameter assignment"""
    X, Y = np.meshgrid(np.linspace(0, 1, 50), np.linspace(0, 1, 50))
    Z = f0 + (1-f0)*X - f0*Y + (2*f0-1)*X*Y
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel('Activator')
    ax.set_ylabel('Repressor')
    ax.set_zlabel('Expression')
    
    # covnert to HTML
    html = mpld3.fig_to_html(fig)
    plt.close(fig)
    return html