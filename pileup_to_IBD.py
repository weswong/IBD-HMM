import sys


def convert_interpreted_pilup_to_ibdtab(file):
    '''analysizes an interpreted bam file and converts to format for steve schaffner's ibd hmm'''
    
    file_out = ('.').join(file.split('.')[0:2]) + '.ibdhmm.txt'
    fout = open(file_out, 'w')
    fin = open(file)
    data= [line.strip().split() for line in fin.readlines()]
    
    ibd_data = [['chrom', 'pos', 'ref', 'a1', 'a2']]
    
    for line in data[1:-1]:
        chr = line[0].split(':')[0].split('_')[1]
        pos = line[0].split(':')[1]
        ref = line[1]
        
        alleles = []
        for i, element in enumerate(line[0:6]):
            if element != '0' and i not in [0,1]:
                alleles.append(data[0][i])
        
        if len(alleles) > 1:
            a1 = alleles[0]
            a2 = alleles[1]
        elif len(alleles) == 0:
            a1, a2 = '.' ,'.'
        else:
            a1, a2 = alleles[0], alleles[0]
        
        new = [str(int(chr)), str(int(pos)), ref, a1, a2]
        ibd_data.append(new)
    for line in ibd_data:
        fout.write(('\t').join(line) + '\n')
        
if __name__ == '__main__':
    file = sys.argv[1]
    convert_interpreted_pilup_to_ibdtab(file)
