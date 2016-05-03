def read(filename):
    import gzip
    import cPickle as pickle
    f = gzip.GzipFile(filename + '.gz', 'rb') 
    obj = pickle.load(f)
    f.close()
    return obj

def write(obj, filename):
    
    import gzip
    import cPickle as pickle
    f = gzip.GzipFile(filename + '.gz', 'wb')
    pickle.dump(obj, f)
    f.close()
    