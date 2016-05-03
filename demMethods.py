def read(filename):
    import bz2
    import cPickle as pickle
    f = bz2.BZ2File(filename + '.bz2', 'rb') 
    obj = pickle.load(f)
    f.close()
    return obj

def write(obj, filename):
    
    import bz2
    import cPickle as pickle
    f = bz2.BZ2File(filename + '.bz2', 'wb')
    s = pickle.dumps(obj)
    for i in range(0, len(s), 2**32):
        f.write(bytes(s[i:i+2**32]))
    f.close()
    