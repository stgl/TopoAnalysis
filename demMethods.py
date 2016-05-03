def read(filename):
    import bz2
    import cPickle as pickle
    f = bz2.BZ2File(filename + '.bz2', 'r') 
    obj = pickle.load(f)
    f.close()
    return obj

def write(obj, filename):
    
    import bz2
    import cPickle as pickle
    f = bz2.BZ2File(filename + '.bz2', 'w')
    pickle.dump(obj, f)
    f.close()
    