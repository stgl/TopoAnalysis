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
    s = pickle.dumps(obj)
    for i in range(0, len(s), 2**32):
        f.write(bytes(s[i:i+2**32]))
    f.close()
    