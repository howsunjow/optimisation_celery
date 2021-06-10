import numpy as np
import pandas as pd
from scipy.optimize import LinearConstraint, Bounds, minimize, curve_fit
from response_functions import *
from collections import namedtuple
from celery import Celery
import celeryconfig


app = Celery('optimiser', broker='pyamqp://guest:guest@157.230.177.180//', backend='rpc://')
app.config_from_object(celeryconfig)
print(app.conf['accept_content'])

class OptimizationParameters(object):
    __slots__ = ['model_parameters', 'initial_adstock', 'costs', 'constraints']
    def __init__(self, model_parameters, initial_adstock, costs, constraints):
        self.model_parameters = model_parameters
        self.initial_adstock = initial_adstock
        self.costs = costs
        self.constraints = constraints

class OptimizationArguments(object):
    __slots__ = ['drivers', 'driver_arguments']
    def __init__(self, drivers, driver_arguments):
        self.drivers = drivers
        self.driver_arguments = driver_arguments

        
FUNCTIONAL_FORMS = {ResponseFunctionType.ADBUDG: 
                        ResponseFunction(adbudg_fun, adbudg_fun_der, adbudg_fun_der2),
                    ResponseFunctionType.LINEAR:
                        ResponseFunction(linear_fun, linear_fun_der, linear_fun_der2),
                    ResponseFunctionType.LOGGED:
                        ResponseFunction(logged_fun, logged_fun_der, logged_fun_der2),
                    ResponseFunctionType.EXPONENTIAL:
                        ResponseFunction(exponential_fun, exponential_fun_der, exponential_fun_der2),
                    ResponseFunctionType.POWER:
                        ResponseFunction(power_fun, power_fun_der, power_fun_der2),
                    ResponseFunctionType.ARCTAN:
                        ResponseFunction(arctan_fun, arctan_fun_der, arctan_fun_der2),
                    ResponseFunctionType.RECIPROCAL:
                        ResponseFunction(reciprocal_fun, reciprocal_fun_der, reciprocal_fun_der2)}

def accumulate_with_initial(inputs, func, initial):
    itr = iter(inputs)
    prev = func(initial, next(itr))
    yield prev
    for cur in itr:
        prev = func(prev, cur)
        yield prev

def create_budget_constraint(num_periods, abudget, num_series=1):
    """This function creates a linear constraint given a number of periods,
        a budget and a number of series. The constraint is such that the spend
        in all periods over all series sums up to the budget

    Args:
        num_periods : int
            The number of distinct periods for which there is spend
        abudget: float
            The budgeted amount that can be spent
        num_series: int
            The number of types of spend

    Returns:
        LinearConstraint: A linear constraint on the spend in all periods
        over all series such that the spend sums up to the budget
    """
    return LinearConstraint([np.ones(num_periods * num_series)],
                            [abudget], [abudget])


def create_spend_bounds(num_periods, abudget, num_series=1):
    """This function creates a bound given a number of periods,
        a budget and a number of series. The bound is such that
        the spend in each period is >= 0 and <= the budget

    Args:
        num_periods : int
            The number of distinct periods for which there is spend
        abudget: float
            The budgeted spend amount
        num_series: int
            The number of types of spend

    Returns:
        Bounds: A Scipy style bound which encodes the spend for each period
        over all series such that each period's spend is bounded between 0 and
        the budgeted spend
    """
    return Bounds(np.zeros(num_periods * num_series),
                  np.ones(num_periods * num_series) * abudget)


def create_plan_zero_period_spend(num_periods, num_series=1):
    return np.zeros(num_periods * num_series)


def get_period_bounds(driver_opt_param, abudget):
    constraints = driver_opt_param.constraints        
    if constraints['period_constraint_type'].lower() == 'units':
        b =  [(0, min(u,abudget)) for u in map(lambda x: constraints['period_maximum']*x , driver_opt_param.costs)]
    else:
        b = [(0, min(abudget,constraints['period_maximum']))] * len(driver_opt_param.costs)     
    r = [(a[0],None if a[1] == np.inf else a[1]) for a in b]
    return(r)

def calculate_total_uplift(plan_spends, args):
    """This calculates the sum of the response functions
       (assumed to be adbug for now) over each plan period,
       over all series taking into account carryover
    Args:
        plan_spend: list
            The list spends for each period and series (spend type)
        driver_parameters: list
            A list of driver parameters

    Returns:
        sum of response functions for all spend periods, across all series
    """
    asum = 0.0
    opt_params = args.driver_arguments
    opt_drivers = args.drivers
    num_series = len(opt_drivers)
    num_periods = len(plan_spends)//num_series
    aseries_idx = 0
    for adriver in opt_drivers:
        # Get the model parameters e.g. coefficient, carryover, etc.
        model_parameters = opt_params.model_parameters[adriver]
        # Get the initial adstocked units at the start of the plan
        initial_adstock = opt_params.initial_adstock[adriver]
        # Get the plan costs
        plan_costs = opt_params.costs[adriver][0:num_periods]
        # From the model parameters get the carryover
        carryover = model_parameters.carryover 
        # Get the offset for the driver spends        
        offset = aseries_idx * num_periods
        # Get the response function
        response_func = FUNCTIONAL_FORMS[model_parameters.response].response
        adstock_units = accumulate_with_initial(
                        map(lambda x,y: x/y , plan_spends[offset:offset+num_periods], plan_costs),
                        lambda x,y: x*carryover + y,
                        initial = initial_adstock) 
        asum += sum(map(lambda x: response_func(x, model_parameters),  adstock_units))
        aseries_idx += 1         
    # Return the negative total uplift so we can minimize   
    return -asum

def calculate_total_uplift_deriv(plan_spends, args):
    """This calculates partial derivative of piecewise_series_adbudg_neg
       wrt to each spend period across all series
    Args:
        plan_spends: list
            The list spends for each period and series (spend type)
        driver_parameters: list
            A list of driver parameters

    Returns:
        a list of partial derivatives for each spend period across all series
    """
    derx = np.zeros(len(plan_spends))
    opt_params = args.driver_arguments
    opt_drivers = args.drivers
    num_series = len(opt_drivers)
    num_periods = len(plan_spends)//num_series
    aseries_idx = 0
    for adriver in opt_drivers:
        # Get the model parameters e.g. coefficient, carryover, etc.
        model_parameters = opt_params.model_parameters[adriver]
        # Get the initial adstocked units at the start of the plan
        initial_adstock = opt_params.initial_adstock[adriver]
        # Get the plan costs
        plan_costs = opt_params.costs[adriver][0:num_periods]
        # From the model parameters get the carryover
        carryover = model_parameters.carryover 
        # Get the offset for the driver spends        
        offset = aseries_idx * num_periods
        
        # Get the response function derivative
        func_deriv = FUNCTIONAL_FORMS[model_parameters.response].derivative
        adstock_units = accumulate_with_initial(
                        map(lambda x,y: x/y , plan_spends[offset:offset+num_periods], plan_costs),
                        lambda x,y: x*carryover + y,
                        initial = initial_adstock)
        adriver_derx = list(map(lambda x,y: func_deriv(x, model_parameters)/y, 
                            adstock_units, plan_costs))
        for i in range(0, num_periods):
            for j in range(i+1, num_periods):
                adriver_derx[i] += adriver_derx[j] * (carryover ** (j-i)) / plan_costs[i]
        for i in range(offset, offset+num_periods):
            derx[i] = -adriver_derx[i-offset]
        aseries_idx += 1                       
    return derx

@app.task
def optimize_budget_spend(abudget,  opt_params, drivers, num_plan_periods = 10):
    num_drivers = len(drivers)
    start_point = create_plan_zero_period_spend(num_plan_periods, num_drivers)
    constraint_list = [create_budget_constraint(num_plan_periods, abudget, num_drivers)]
    opt_args = OptimizationArguments(drivers, opt_params)
    res = minimize(calculate_total_uplift, start_point,
                args=opt_args, method='SLSQP', constraints=constraint_list,
                bounds=[(0, None) for i in range(len(start_point))],
                options={'disp': True, 'maxiter': 1000})
    return(res)
