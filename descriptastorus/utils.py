def raw_to_libsvm(ofile, store, full_first_header=False):
    """Write a raw store to libsvm format
    
     :param ofile: file like object to write to
     :param store: the raw store
     :param full_first_header: Write the full first line regardless of sparseness
    """
    
    for idx,data in enumerate(store):
        line = []
        for label_idx, (c, v) in enumerate(zip(store.colnames, data)):
            if full_first_header and idx == 0:
                line.append("%s:%f"%(label_idx, v))
            elif v:
                line.append("%s:%f"%(label_idx, v))
        ofile.write(" ".join(line))
        ofile.write("\n")

        
