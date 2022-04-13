def parse_genome(genome_file):
    """
    """
    genome_dict = {}

    header = ""
    seq = ""

    with open(genome_file) as f:  # read genome file
        for line in f:  # parse genome file
            if line.startswith(">"):  # parse header
            #* if there is a header already, we store the current sequence
            #* to this header.
                if header:
                    genome_dict[header] = seq
                    #* then we flush the sequence
                    seq = ""
                #* and parse the new header
                header = line.strip().split(" ")[0][1:]  # remove '>' and ' '
            #elif line != "\n":  # parse sequence and skips newlines
            else:
                #* if no header is detected, we append the new line
                #* to the current sequence
                seq += line.strip()
                #genome_dict[name] = line.strip()  # append sequence to entry
        #* after the last line, we have to store
        #* the last sequence as well. Since no new line with ">" occurs, we
        #* do this manually
        genome_dict[header] = seq
    return genome_dict