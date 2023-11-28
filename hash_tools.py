from hashlib import sha256
import numpy
import scipy

def hash(array):
    if isinstance(array, list):
        if isinstance(array[0], list):
            concat = []
            for sublist in array:
                concat += sublist
                #mark the end of the sublist with a distinctive feature
                concat.append(len(sublist))
            to_hash = numpy.array(concat)
        else:
            to_hash = numpy.array(array)
    elif scipy.sparse.issparse(array):
        to_hash = array.todense()
    elif isinstance(array, numpy.ndarray):
        to_hash = array
    else:
        raise TypeError("Unsupported type", type(array))

    return sha256(to_hash.tobytes()).hexdigest()
