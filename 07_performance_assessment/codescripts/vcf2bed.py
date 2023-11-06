import re

def position_to_bed(chromo,pos,ref,alt):
    # snv
    # Start cooridinate BED = start coordinate VCF - 1
    # End cooridinate BED = start coordinate VCF 

    if len(ref) == 1 and len(alt) == 1:
        StartPos = int(pos) -1
        EndPos = int(pos)
    
    # deletions
    # Start cooridinate BED = start coordinate VCF - 1
    # End cooridinate BED = start coordinate VCF + (reference length - alternate length)

    elif len(ref) > 1 and len(alt) == 1:
        StartPos = int(pos) - 1
        EndPos = int(pos) + (len(ref) - 1)
        
    #insertions
    # For insertions:
    # Start cooridinate BED = start coordinate VCF - 1
    # End cooridinate BED = start coordinate VCF + (alternate length - reference length)

    else:
        StartPos = int(pos) - 1
        EndPos = int(pos) + (len(alt) - 1)

    return chromo,StartPos,EndPos

def padding_region(chromo,pos1,pos2,padding):
    StartPos1 = pos1 - padding
    EndPos1 = pos1
    StartPos2 = pos2
    EndPos2 = pos2 + padding
    return chromo,StartPos1,EndPos1,StartPos2,EndPos2
