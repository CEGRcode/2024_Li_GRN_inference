import copy
from scipy.integrate import solve_ivp
import numpy as np
import sympy as sp
from scipy.optimize import minimize
import random
from scipy import stats
from scipy.optimize import minimize_scalar
import re
import numexpr as ne
from scipy.stats import norm

def fit_n_from_distribution_list_return_n(func_str,
                                         dist_list,
                                         weights=None,
                                         grid_size_side=100,   # grid_side x grid_side -> grid_side**2 points
                                         bins=100,
                                         n_bounds=(0.01, 1000.0),
                                         z_range=(0.0, 1.0),
                                         tol=1e-3,
                                         verbose=False):
    """
    Fit scalar n so that the histogram of Z = F(X,Y;n) on an XY grid matches
    a target mixture of normals specified by dist_list = [(mu1,sigma1), ...].

    Returns:
        n_hat (float) -- fitted n
    """
    # ---- helpers (nested so this is a single copy-paste function) ----
    def rewrite_hill_safe(s):
        # replace X**n/(X**n + C**n) -> 1/(1 + (C/X)**n) and same for Y
        pat_X = re.compile(r'(?P<L>X)\s*\*\*\s*n\s*/\s*\(\s*X\s*\*\*\s*n\s*\+\s*(?P<C>[0-9eE\.\+\-]+)\s*\*\*\s*n\s*\)')
        s = pat_X.sub(r'1.0/(1.0 + (\g<C>/\1)**n)', s)
        pat_X2 = re.compile(r'\(\s*X\s*\*\*\s*n\s*\)\s*/\s*\(\s*X\s*\*\*\s*n\s*\+\s*(?P<C>[0-9eE\.\+\-]+)\s*\*\*\s*n\s*\)')
        s = pat_X2.sub(r'1.0/(1.0 + (\g<C>/X)**n)', s)
        pat_Y = re.compile(r'(?P<L>Y)\s*\*\*\s*n\s*/\s*\(\s*Y\s*\*\*\s*n\s*\+\s*(?P<C>[0-9eE\.\+\-]+)\s*\*\*\s*n\s*\)')
        s = pat_Y.sub(r'1.0/(1.0 + (\g<C>/\1)**n)', s)
        pat_Y2 = re.compile(r'\(\s*Y\s*\*\*\s*n\s*\)\s*/\s*\(\s*Y\s*\*\*\s*n\s*\+\s*(?P<C>[0-9eE\.\+\-]+)\s*\*\*\s*n\s*\)')
        s = pat_Y2.sub(r'1.0/(1.0 + (\g<C>/Y)**n)', s)
        return s

    def evaluate_string_on_grid_safe(s, n_val, X_arr, Y_arr):
        loc = {'X': X_arr, 'Y': Y_arr, 'n': n_val}
        try:
            Z = ne.evaluate(s, local_dict=loc)
            return np.asarray(Z, dtype=float).ravel()
        except Exception:
            # fallback eval with a safe small environment and a stable hill if needed
            def hill(u, n, t):
                u = np.asarray(u, dtype=float)
                eps = 1e-12
                ratio = np.where(u > eps, (t / (u + eps))**n, np.inf)
                return np.where(u > eps, 1.0 / (1.0 + ratio), 0.0)
            safe_locals = {'X': X_arr, 'Y': Y_arr, 'n': n_val, 'hill': hill, 'np': np}
            Z = eval(s, {'__builtins__': None}, safe_locals)
            return np.asarray(Z, dtype=float).ravel()

    def hist_density(samples, bins, z_min, z_max):
        counts, edges = np.histogram(samples, bins=bins, range=(z_min, z_max), density=False)
        widths = np.diff(edges)
        density = counts / (np.sum(counts * widths) + 1e-16)
        return 0.5*(edges[:-1]+edges[1:]), density

    def make_target_pdf_from_normals(dist_list, weights, z_min, z_max, bins):
        k = len(dist_list)
        if weights is None:
            w = np.ones(k) / k
        else:
            w = np.asarray(weights, dtype=float)
            w = w / np.sum(w)
        edges = np.linspace(z_min, z_max, bins+1)
        centers = 0.5*(edges[:-1] + edges[1:])
        pdf = np.zeros_like(centers, dtype=float)
        for (mu, sigma), wi in zip(dist_list, w):
            if sigma <= 0:
                sigma = 1e-6
            pdf += wi * norm.pdf(centers, loc=mu, scale=sigma)
        # renormalize over [z_min,z_max]
        area = np.trapz(pdf, centers)
        if area <= 0:
            raise ValueError("Target PDF has zero area; check dist_list/bins.")
        pdf /= area
        return centers, pdf

    # ---- prepare XY grid ----
    z_min, z_max = z_range
    xs = np.linspace(0.0, 1.0, grid_size_side)
    ys = np.linspace(0.0, 1.0, grid_size_side)
    Xg, Yg = np.meshgrid(xs, ys)
    Xflat = Xg.ravel()
    Yflat = Yg.ravel()

    # ---- sanitize func_str for numeric stability ----
    func_safe = rewrite_hill_safe(func_str)

    # ---- build target pdf (on histogram centers) ----
    z_centers_target, pdf_target = make_target_pdf_from_normals(dist_list, weights, z_min, z_max, bins)

    # ---- objective: L2 between model histogram and target pdf ----
    def objective(n_scalar):
        n = float(n_scalar)
        if n <= 0:
            return 1e9 + (1.0 - n)**2
        try:
            Z = evaluate_string_on_grid_safe(func_safe, n, Xflat, Yflat)
        except Exception as e:
            if verbose:
                print("eval error at n=", n, ":", e)
            return 1e6 + abs(n)*1e3
        Z = np.clip(Z, z_min, z_max)
        _, pdf_model = hist_density(Z, bins=bins, z_min=z_min, z_max=z_max)
        resid = pdf_model - pdf_target
        return float(np.sum(resid**2))

    # ---- optimize scalar n ----
    res = minimize_scalar(objective, bounds=n_bounds, method='bounded', options={'xatol': tol})
    n_hat = float(res.x)
    return n_hat


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

def Run_Dynamics(GRN_instance_ori, index_i_func, index_j_func, TrainningCount_func, WTTP, TRMAX, return_dict, index_of_diff_genes):
    '''Simulate dynamics'''
    t_ticks = [int(TrainningCount_func/250)*tick for tick in range(0, 250+1)]
    GRN_instance = copy.deepcopy(GRN_instance_ori)
    specific_gene = int(GRN_instance.Name.split('_')[1])
    GRN_instance.SetmRNA(WTTP[str(index_i_func)][1], WTTP[str(index_i_func)][0], TRMAX)
    GRN_instance.Exclude_Non_Diff_Genes(index_of_diff_genes)
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
    sol = [x.item() for x in list(sol.T)]
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