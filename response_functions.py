"""
This module contains variety of response functions (and their derivatives)
"""
import math
from collections import namedtuple
from enum import Enum


class ResponseFunctionType(Enum):
    ADBUDG = 1
    LINEAR = 2
    LOGGED = 3
    EXPONENTIAL = 4
    POWER = 5
    ARCTAN = 6
    RECIPROCAL = 7

Parameter = namedtuple('Parameter',
                       ['response', 'coefficient', 'gamma', 'rho', 'carryover', 'lag'])
Parameter.__new__.__defaults__ = (ResponseFunctionType.ADBUDG, 0.0, math.nan, math.nan, 0.0, 0)

ResponseFunction = namedtuple('ResponseFunction',
                              ['response', 'derivative', 'second_derivative'])


def adbudg_fun(spend: float, adbudg_parameters: Parameter):
    """The Adbudg response function

    Args:
        spend : double
            The value at which to calculate the response
        adbudg_parameters: Parameter
            adbug parameters as a named tuple

    Returns:
        double: The response to the input spend
    """
    beta = adbudg_parameters.coefficient
    gamma = adbudg_parameters.gamma
    rho = adbudg_parameters.rho
    return (beta * spend**gamma)/(spend**gamma + rho**gamma)


def adbudg_fun_der(spend: float, adbudg_parameters: Parameter):
    """The derivative of the Adbudg response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        adbudg_parameters: Parameter
            adbug parameters as a named tuple

    Returns:
        float: The derivative of the response at spend

    """
    beta = adbudg_parameters.coefficient
    gamma = adbudg_parameters.gamma
    rho = adbudg_parameters.rho
    numerator = beta * gamma * rho**gamma * spend**(gamma-1.0) 
    denom = (spend**gamma + rho**gamma)**2.0
    return numerator/denom


def adbudg_fun_der2(spend: float, adbudg_parameters: Parameter):
    """The second derivative of the Adbudg response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        adbudg_parameters: Parameter
            adbug parameters as a named tuple

    Returns:
        float: The second derivative of the response at spend

    """
    beta = adbudg_parameters.coefficient
    gamma = adbudg_parameters.gamma
    rho = adbudg_parameters.rho
    return ((beta * rho**gamma)/((spend**gamma + rho**gamma)**3)) * \
           ((-gamma-1) * spend**(2*(gamma-1)) +
            (rho**gamma * (gamma-1) * spend**(gamma-2)))


def linear_fun(spend: float, linear_parameters: Parameter):
    """The linear response function

    Args:
        spend : double
            The value at which to calculate the response
        linear_parameters: Parameter
            Function parameters

    Returns:
        double: The response to the input spend

    """
    beta = linear_parameters.coefficient
    return beta * spend


def linear_fun_der(_: float, linear_parameters: Parameter):
    """The derivative of the linear response function

    Args:
        spend : double
            The value at which to calculate the response
        linear_parameters: Parameter
            Function parameters


    Returns:
        double: The derivative of the response at spend

    """
    return linear_parameters.coefficient


def linear_fun_der2(_spend: float, _linear_parameters: Parameter):
    """The second derivative of the linear response function

    Args:
        spend : double
            The value at which to calculate the response
        linear_parameters: Parameter
            Function parameters


    Returns:
        double: The derivative of the response at spend

    """
    return 0.0


def logged_fun(spend: float, logged_parameters: Parameter):
    """The logged response function

    Args:
        spend : double
            The value at which to calculate the response
        logged_parameters: Parameter
            Function parameters

    Returns:
        double: The response to the input spend

    """
    return logged_parameters.coefficient * math.log(spend)


def logged_fun_der(spend: float, logged_parameters: Parameter):
    """The derivative of the logged response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        logged_parameters: Parameter
            Function parameters

    Returns:
        double: The derivative of the response at spend

    """
    return logged_parameters.coefficient/spend


def logged_fun_der2(spend: float, logged_parameters: Parameter):
    """The second derivative of the logged response function

    Args:
        spend : double
            The value at which to calculate the second derivative
            of the response
        logged_parameters: Parameter
            Function parameters

    Returns:
        double: The second derivative of the response at spend

    """
    return -logged_parameters.coefficient/(spend**2)


def exponential_fun(spend: float, exponential_parameters: Parameter):
    """The exponential (diminishing returns) response function

    Args:
        spend : double
            The value at which to calculate the response
        exponential_parameters: Parameter
            Function parameters

    Returns:
        double: The response to the input spend

    """
    beta = exponential_parameters.coefficient
    gamma = exponential_parameters.gamma
    return beta * (1 - math.exp(-spend/gamma))


def exponential_fun_der(spend: float, exponential_parameters: Parameter):
    """The derivative of the exponential (diminishing returns) response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        exponential_parameters: Parameter
            Function parameters

    Returns:
        double: The derivative of the response at spend

    """
    beta = exponential_parameters.coefficient
    gamma = exponential_parameters.gamma
    return beta/gamma * math.exp(-spend/gamma)


def exponential_fun_der2(spend: float, exponential_parameters: Parameter):
    """The second derivative of the exponential (diminishing returns)
       response function

    Args:
        spend : double
            The value at which to calculate the second derivative of
            the response
        exponential_parameters: Parameter
            Function parameters

    Returns:
        double: The second derivative of the response at spend

    """
    beta = exponential_parameters.coefficient
    gamma = exponential_parameters.gamma
    return -beta/(gamma**2) * math.exp(-spend/gamma)


def power_fun(spend: float, power_parameters: Parameter):
    """The power response function

    Args:
        spend : double
            The value at which to calculate the response
        power_parameters: Parameter
            Function parameters


    Returns:
        double: The response to the input spend

    """
    beta = power_parameters.coefficient
    gamma = power_parameters.gamma
    return beta * (spend**gamma)


def power_fun_der(spend: float, power_parameters: Parameter):
    """The derivative of the power response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        power_parameters: Parameter
            Function parameters

    Returns:
        double: The derivative of the response at spend

    """
    beta = power_parameters.coefficient
    gamma = power_parameters.gamma
    return beta * gamma * (spend**(gamma-1))


def power_fun_der2(spend: float, power_parameters: Parameter):
    """The second derivative of the power response function

    Args:
        spend : double
            The value at which to calculate the second derivative
            of the response
        power_parameters: Parameter
            Function parameters

    Returns:
        double: The second derivative of the response at spend

    """
    beta = power_parameters.coefficient
    gamma = power_parameters.gamma
    return beta * gamma * ((gamma-1.0) * spend**(gamma-2.0))


def arctan_fun(spend: float, arctan_parameters: Parameter):
    """The arctan response function

    Args:
        spend : double
            The value at which to calculate the response
        arctan_parameters: Parameter
            Function parameters

    Returns:
        double: The response to the input spend

    """
    beta = arctan_parameters.coefficient
    gamma = arctan_parameters.gamma
    return (2.0 * beta/math.pi) * math.atan(spend/gamma)


def arctan_fun_der(spend: float, arctan_parameters: Parameter):
    """The derivative of the arctan response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        arctan_parameters: Parameter
            Function parameters

    Returns:
        double: The derivative of the response at spend

    """
    beta = arctan_parameters.coefficient
    gamma = arctan_parameters.gamma
    return (2.0 * beta * gamma) / (math.pi * (gamma**2 + spend**2))


def arctan_fun_der2(spend: float, arctan_parameters: Parameter):
    """The second derivative of the arctan response function

    Args:
        spend : double
            The value at which to calculate the second derivative
            of the response
        arctan_parameters: Parameter
            Function parameters

    Returns:
        double: The second derivative of the response at spend

    """
    beta = arctan_parameters.coefficient
    gamma = arctan_parameters.gamma
    numerator = (-2.0 * beta * gamma) * (2.0 * spend)
    denom = (math.pi * ((gamma**2 + spend**2)**2))
    return numerator/denom


def reciprocal_fun(spend: float, reciprocal_parameters: Parameter):
    """The reciprocal response function

    Args:
        spend : double
            The value at which to calculate the response
        reciprocal_parameters: Parameter
            Function parameters

    Returns:
        double: The response to the input spend

    """
    beta = reciprocal_parameters.coefficient
    rho = reciprocal_parameters.rho
    return beta/(spend + rho/10.0)


def reciprocal_fun_der(spend: float, reciprocal_parameters: Parameter):
    """The derivative of the reciprocal response function

    Args:
        spend : double
            The value at which to calculate the derivative of the response
        reciprocal_parameters: Parameter
            Function parameters

    Returns:
        double: The derivative of the response at spend

    """
    beta = reciprocal_parameters.coefficient
    rho = reciprocal_parameters.rho
    return -beta/((spend + rho/10.0)**2)


def reciprocal_fun_der2(spend: float, reciprocal_parameters: Parameter):
    """The second derivative of the reciprocal response function

    Args:
        spend : double
            The value at which to calculate the second derivative of the
            response
        reciprocal_parameters: Parameter
            Function parameters

    Returns:
        double: The second derivative of the response at spend

    """
    beta = reciprocal_parameters.coefficient
    rho = reciprocal_parameters.rho
    return 2.0 * beta/((spend + rho/10)**3)

ResponseFunction.__new__.__defaults__ = (adbudg_fun, adbudg_fun_der, adbudg_fun_der2)
