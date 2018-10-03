import numpy as np
import matplotlib.pyplot as plt
import inspect
from functools import wraps
import csv
import sys

SAMPLE_CSV = "d11B-pH_Error_estimation.csv"
SAMPLE_OUT_CSV = "d11B-pH_Error_estimation_out.csv"

def showparam(f):
    @wraps(f)
    def wrapper(*args, **kwds):
        arg_names, varargs, keywords, defaults = inspect.getargspec(f)
        print("Parameters for function {n}:".format(n=f.__name__))
        for arg, arg_val in zip(arg_names, args):
            print("    {a}: {v}".format(a=arg, v=arg_val))
        return f(*args, **kwds)
    return wrapper

ALPHA = 1.0272
N_SAMPLE = 1000000
FIELDNAMES = ['Sample ID', 'Salinity', 'SD of Salinity', 'Temperature', 'SD of Temperature', 'pKb', 'SD of pKb', 'd11Bsw', 'SD of d11Bsw', 'd11Bc', 'SD of d11Bc', 'pHcf', 'SD of pHcf']
OUTPUT_FIELDS = ['pKb', 'SD of pKb', 'pHcf', 'SD of pHcf']
MENDATORY_FIELDS = ['Salinity', 'SD of Salinity', 'Temperature', 'SD of Temperature', 'd11Bsw', 'SD of d11Bsw', 'd11Bc']
CALC_FIELDS = MENDATORY_FIELDS + ['SD of d11Bc']


def print_stats(var, var_name="unnamed_variable"):
    print("{v:12}: avg{m:12}, stdev={s:12}".format(v=var_name, m=var.mean(), s=var.std()))

def plot_hist(var, var_name="unnamed_variable", bin_count=30):
    count, bins, ignored = plt.hist(var, bin_count, label=var_name)
    plt.plot(bins)
    plt.show()
    
def present_var(var, var_name="unnamed_variable"):
    """Presents one variable only
    """
    print_stats(var, var_name)
    plot_hist(var, var_name)
    
def present_vars(vars, var_names, bin_count=30):
    """Presents multiple variables on a single chart
    """
    for var, var_name in zip(vars, var_names):
        print_stats(var, var_name)
        plt.hist(var, bin_count, range=(min(var), max(var)), alpha=0.5, label=var_name)
    plt.legend(loc='upper right')
    plt.show()

def demo():
    """Demo with formula: y = x^2 + 2*x + 1
    Where x is a random variable that follows normal distribution with mean=0 & sigma=1
    """
    mu, sigma = 0, 1
    print("N={n}, mu={m}, sigma={s}".format(n=N_SAMPLE, m=mu, s=sigma))
    x = sampling(mu, sigma, N_SAMPLE)   # N-sized array of random variable x
    y = x**2 + 2*x + 1                          # N-sized array of deduced variable y
    present_vars((x, y), ("x", "y"))

def pH_formula(pKb, d11Bsw, d11Bcarbonate, alpha):
    return pKb - np.log10((d11Bsw - d11Bcarbonate) / (alpha * d11Bcarbonate - d11Bsw + 1000 * (alpha - 1)))

@showparam
def calc_pH_by_d11Bcarbonate(mu_d11Bcarbonate, sigma_d11Bcarbonate, N=N_SAMPLE):
    # random variable
    d11Bcarbonate = sampling(mu_d11Bcarbonate, sigma_d11Bcarbonate, N)
    # constants
    pKb = 8.6152
    d11Bsw = 39.61
    alpha = 1.0272
    print("pKb = {k}, d11Bsw = {s}, alpha = {a}".format(k=pKb, s=d11Bsw, a=alpha))
    # formula
    pH = pH_formula(pKb, d11Bsw, d11Bcarbonate, alpha)
    present_vars((d11Bcarbonate, pH), ("d11Bcarbonate", "pH"))

def gen_pKb_by_S_and_T(mu_pKb_S, sigma_pKb_S, mu_pKb_T, sigma_pKb_T, N=N_SAMPLE):
    S = sampling(mu_pKb_S, sigma_pKb_S, N)
    T = sampling(mu_pKb_T+273.15, sigma_pKb_T, N)
#     k1, k2, k3, k4, k5, k6 = 8966.9, 2890.53, 77.942, 1.728, 0.0996, 148.0248, 137.1942, 1.62142, 24.4344, 25.085, 0.2474, 0.053105
    return -np.log10(np.exp((-8966.9-2890.53*(S**0.5)-77.942*S+1.728*(S**1.5)-0.0996*(S**2))/T+148.0248+137.1942*(S**0.5)+1.62142*S-(24.4344+25.085*(S**0.5)+0.2474*S)*np.log(T)+0.053105*(S**0.5)*T))

def sampling(mu, std, N):
    if std < 0:
        raise ValueError
    elif std == 0:
        return np.array([mu]*N)
    else:
        return np.random.normal(mu, std, N)

def gen_valid_sample(s, s_std, t, t_std, sw, sw_std, c, c_std, N=N_SAMPLE):
    d11Bcarbonate = sampling(c, c_std, N)
    pKb = gen_pKb_by_S_and_T(s, s_std, t, t_std, N)
    d11Bsw = sampling(sw, sw_std, N)
    alpha = 1.0272
    pH = pH_formula(pKb, d11Bsw, d11Bcarbonate, alpha)
    pH_valid_mask = ~np.isnan(pH)
    valid_ratio = 1.0*len(pH[pH_valid_mask])/N
    return d11Bcarbonate[pH_valid_mask], pKb[pH_valid_mask], d11Bsw[pH_valid_mask], pH[pH_valid_mask], valid_ratio

@showparam
def calc_pH_by_d11Bcarbonate_with_random_parameters(mu_d11Bcarbonate, sigma_d11Bcarbonate, mu_pKb_S, sigma_pKb_S, mu_pKb_T, sigma_pKb_T, mu_d11Bsw, sigma_d11Bsw, N=N_SAMPLE):
    d11Bcarbonate_valid, pKb_valid, d11Bsw_valid, pH_valid, valid_ratio = gen_valid_sample(mu_pKb_S, sigma_pKb_S, mu_pKb_T, sigma_pKb_T, mu_d11Bsw, sigma_d11Bsw, mu_d11Bcarbonate, sigma_d11Bcarbonate, N)
    print("valid sample count: {v} ({p}%)".format(v=len(pH_valid), p=100.0*valid_ratio))
    present_vars((d11Bcarbonate_valid, pKb_valid, d11Bsw_valid, pH_valid), ("d11Bcarbonate", "pKb", "d11Bsw", "pH"))
      


DEFAULT_SD_OF_D11BC = 0.82
      
def process_csv(filename, outname):
    print(filename)
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        if reader.fieldnames != FIELDNAMES:
            print("error: illegal fieldnames: {n}".format(n=reader.fieldnames))
        data = [row for row in reader]
    with open(outname, 'wb') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in data:
#             print row
            if all([row[k] for k in OUTPUT_FIELDS]):
                print("skipped Sample ID {s} which seems already done".format(s=row["Sample ID"]))
            else:
                try:
                    assert all([row[k] for k in MENDATORY_FIELDS])
                    if not row['SD of d11Bc']:
                        row['SD of d11Bc'] = DEFAULT_SD_OF_D11BC
                    for field in CALC_FIELDS:
    #                     print field, row[field]
                        row[field] = float(row[field])
                    d11Bcarbonate_valid, pKb_valid, d11Bsw_valid, pH_valid, valid_ratio = gen_valid_sample(*[row[k] for k in CALC_FIELDS])
                    row['pKb'] = pKb_valid.mean()
                    row['SD of pKb'] = pKb_valid.std()
                    row['pHcf'] = pH_valid.mean()
                    row['SD of pHcf'] = pH_valid.std()
                    print("Sample ID {s:16}: {k:8.5f}, {t:8.5f}, {h:8.5f}, {d:8.5f} ratio={r:7.3f}%".format(
                        s=row["Sample ID"], 
                        k=row['pKb'], 
                        t=row['SD of pKb'], 
                        h=row['pHcf'], 
                        d=row['SD of pHcf'],
                        r=100.0*valid_ratio))
                except Exception as e:
                    print(e)
                    print("warning: skipped Sample ID {s}".format(s=row["Sample ID"]))
            writer.writerow(row)

if __name__ == "__main__":
    # Run demo
    # Replace with actual formula for real use
#     demo()
    
#     calc_pH_by_d11Bcarbonate(17.39, 0.19)

#     calc_pH_by_d11Bcarbonate_with_random_parameters(12.35, 0.66/2, 32, 0.2, 25, 0.1, 39.61, 0.4/2)

    if len(sys.argv) == 1:
        csvname = SAMPLE_CSV
        outname = SAMPLE_OUT_CSV
    else:
        csvname = sys.argv[1]
        outname = csvname
    process_csv(csvname, outname)