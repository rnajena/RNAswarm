/*************************************************************************
* annotation table to fasta
**************************************************************************/
process annotationTableToFasta {
    label: "RNAswarm_small"

    input:
    tuple
}


/*************************************************************************
* predict RNA-RNA structures
**************************************************************************/
process runRNAcofold {
    label: "viennarna"

    input:
    tuple 
}

// import os

// # Parse pretty table, which is a tsv file
// pretty_table = '/Users/gabriellovate/Research/Presentations/20240220-ITN_update/data/pretty_table.tsv'
// with open(pretty_table, 'r') as f:
//     for line in f:
//         # skip the header
//         if line.startswith('Unnamed: 0'):
//             continue
//         # split the line into fields
//         id, segment01, start01, end01, segment02, start02, end02, area, mean_by_area, overlaps_with, mean_readcount, seq01, seq02, seq01_extended, seq02_extended, smiley, position_seq01, position_seq02, energy, smiley_extended, position_seq01_extended, position_seq02_extended, energy_extended = line.strip().split('\t')
//         # Check if overlaps_with is not empty, if not, skip
//         if overlaps_with:
//             continue
//         start_seq01, end_seq01 = position_seq01.split(',')
//         start_seq02, end_seq02 = position_seq02.split(',')
//         start_seq01_extended, end_seq01_extended = position_seq01_extended.split(',')
//         start_seq02_extended, end_seq02_extended = position_seq02_extended.split(',')
//         # Convert to 0-based index exclusive
//         start_seq01, end_seq01 = int(start_seq01) - 1, int(end_seq01)
//         start_seq02, end_seq02 = int(start_seq02) - 1, int(end_seq02)
//         start_seq01_extended, end_seq01_extended = int(start_seq01_extended) - 1, int(end_seq01_extended)
//         start_seq02_extended, end_seq02_extended = int(start_seq02_extended) - 1, int(end_seq02_extended)
//         # Get the sequences
//         seq01 = seq01[start_seq01:end_seq01]
//         seq02 = seq02[start_seq02:end_seq02]
//         seq01_extended = seq01_extended[start_seq01_extended:end_seq01_extended]
//         seq02_extended = seq02_extended[start_seq02_extended:end_seq02_extended]
//         # Run RNAplot to generate svg files of the interactions
//         os.system(f'echo ">{id}\n{seq01}&{seq02}\n{smiley}" | RNAplot -o svg')
//         os.system(f'echo ">{id}extended\n{seq01_extended}&{seq02_extended}\n{smiley_extended}" | RNAplot -o svg')
//         # Run RNAplot to generate ps files of the interactions
//         os.system(f'echo ">{id}\n{seq01}&{seq02}\n{smiley}" | RNAplot -o ps')
//         os.system(f'echo ">{id}extended\n{seq01_extended}&{seq02_extended}\n{smiley_extended}" | RNAplot -o ps')