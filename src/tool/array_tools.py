from typing import Literal
import numpy as np

def tuple_from_array(arr: np.ndarray | tuple, expected_lengths: int | list[int], desired_len) -> tuple[tuple, int, str]:
    """
    Convert a numpy array or tuple to a specified format.

    Parameters
    ----------
    arr : Union[np.ndarray, Tuple]
        Input array or tuple.
    expected_lengths : Union[int, List[int]]
        Expected length(s) of the input array. If a single integer is provided,
        it will be wrapped in a list.
    desired_len : int
        Desired length of the output tuple.

    Returns
    -------
    Tuple[Tuple, int, str]
        A tuple containing:
            - The converted values as a tuple.
            - The original input length.
            - The type of input ("array" or "tuple").

    Raises
    ------
    RuntimeError
        If the input does not match the expected lengths.

    """
    expected_lengths = expected_lengths if isinstance(expected_lengths, list) else [expected_lengths]

    if isinstance(arr, tuple):
        if len(arr) not in expected_lengths:
            raise RuntimeError("Error in array_tools.tuple_from_array(): " +
                               f"expected length of array to be one of {expected_lengths}, but is {len(arr)}!\n\t{arr=}")
        input_type = "tuple"
        input_len = len(arr)
        tvals = arr

    else:
        if arr.shape[0] == 1:
            arr = arr.transpose()
        if arr.shape[0] not in expected_lengths:
            raise RuntimeError("Error in array_tools.tuple_from_array(): " +
                               f"expected length of array to be one of {expected_lengths}, but is {arr.shape[0]}!\n\t{arr=}")
        if arr.ndim > 1 and arr.shape[1] > 1:
            raise RuntimeError("Error in array_tools.tuple_from_array(): " +
                               f"Can only currently accept single-length arrays (shape Nx1)!\n\t{arr=}")

        input_type = "array"
        input_len = arr.shape[0]
        tvals = tuple(arr.tolist())

    if len(tvals) == desired_len:
        return tvals, input_len, input_type
    else:
        return_vals = [0]*desired_len
        for i in range(min(input_len, desired_len)):
            return_vals[i] = tvals[i]
        return tuple(return_vals), input_len, input_type

def retval_from_tuple(retval: np.ndarray | tuple, return_len: int, return_type: Literal["array"] | Literal["tuple"]) -> np.ndarray | tuple:
    """
    Convert a numpy array or tuple to the specified output type and length.

    Parameters
    ----------
    retval : Union[np.ndarray, Tuple]
        Input array or tuple.
    return_len : int
        Desired length of the output array/tuple.
    return_type : Literal["array"] | Literal["tuple"]
        The desired output type ("array" or "tuple").

    Returns
    -------
    Union[np.ndarray, Tuple]
        The converted values in the specified format.

    Raises
    ------
    ValueError
        If return_type is not one of "array" or "tuple".

    """
    if isinstance(retval, tuple):
        input_len = len(retval)
    else:
        if retval.shape[0] == 1:
            retval = retval.transpose()
        input_len = retval.shape[0]

    tvals, _, _ = tuple_from_array(retval, input_len, return_len)

    if return_type == "tuple":
        return tvals
    elif return_type == "array":
        return np.array(tvals)
    else:
        raise ValueError("Error in array_tools.retval_from_tuple(): " +
                         f"return_type must be one of \"array\" or \"tuple\", but is {return_type}!")
