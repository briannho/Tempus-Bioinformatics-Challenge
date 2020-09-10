#!/usr/bin/env python3

# import libraries
import argparse
from argparse import RawTextHelpFormatter
from pysam import VariantFile
import requests
import json

# returns parameters
def get_args():

    # define parameters
    argparser = argparse.ArgumentParser(description = 'This is a prototype VCF annotation tool. \
                                        \nexample: python3 tempus_bioi_challenge.py -i /path/to/in.vcf -o /path/to/out.tsv', \
                                        formatter_class=RawTextHelpFormatter)
    argparser.add_argument("-i", "--input", required = True, help = "vcf file")
    argparser.add_argument("-o", "--output", required = True, help = "output tsv file")
    args = argparser.parse_args()
    
    return args

# fetch variants from vcf file
def get_variants(vcf):
    return [var for var in vcf.fetch()]
    
# fetch sample names from vcf file
def get_samples(vcf):
    return list(vcf.header.samples)

# returns variant type(s) for each variant
def variant_types(variants):
    var_types = []
    
    for var in variants:
        var_type = var.info['TYPE']
        
        # one type associated with variant
        if len(var_type) == 1:
            var_types.append(var_type[0])
            
        # multiple types associated with variant because there are multiple alleles
        elif len(var_type) > 1:
            var_types.append(','.join(var.info['TYPE']))
            
    print('Variant types annotated.')
    return var_types

# returns effects caused by variants from ExAC API
def effects(variants):
    effects = []
    
    for var in variants:
        chrom = str(var.chrom)
        pos = str(var.pos)
        ref = var.ref
        alts = var.alts

        # one allele associated with variant
        if len(alts) == 1:

            # ExAC API request
            var_id = f'{chrom}-{pos}-{ref}-{alts[0]}'
            request = f'http://exac.hms.harvard.edu/rest/variant/consequences/{var_id}'
            response = requests.get(request)
            text = json.loads(response.text)
            
            # No recorded effects associated with variant
            if text is None or text == {}:
                effects.append('.')
            
            # there is at least one effect associated  
            else:
                effect = list(text.keys())
                
                # only one effect is associated
                if len(effect) == 1:
                    effects.append(str(effect[0]))
                
                # more than one effects are associated, so the most deleterious effect is chosen
                elif len(effect) > 1:
                    effects.append(rank_effects(effect))

        # more than one alleles associated with variant
        if len(alts) > 1:
            effects_temp = []
            
            # iterate through each allele
            for alt in alts:

                # ExAC API request
                var_id = f'{chrom}-{pos}-{ref}-{alt}'
                request = f'http://exac.hms.harvard.edu/rest/variant/consequences/{var_id}'
                response = requests.get(request)
                text = json.loads(response.text)
                
                # No recorded effects associated with variant
                if text is None or text == {}:
                    effects_temp.append('.')
                
                # there is at least one effect associated
                else:
                    effect = list(text.keys())
                    
                    # only one effect is associated
                    if len(effect) == 1:
                        effects_temp.append(str(effect[0])) 
                    
                    # more than one effects are associated, so the most deleterious effect is chosen
                    elif len(effect) > 1:
                        effects_temp.append(rank_effects(effect))

            # concatenates elements of effects_temp into one string and stores it in effects
            effects.append(','.join(effects_temp))
    
    print('Variant effects annotated.')
    return effects

# returns most deleterious effect if variant has two or more effects
def rank_effects(effects):
    if 'stop_gained' in effects:
        return 'stop_gained'
    elif 'stop_lost' in effects:
        return 'stop_lost'
    elif 'missense_variant' in effects:
        return 'missense_variant'
    elif 'initiator_codon_variant' in effects:
        return 'initiator_codon_variant'
    elif 'splice_acceptor_variant' in effects:
        return 'splice_acceptor_variant'
    elif 'splice_donor_variant' in effects:
        return 'splice_donor_variant'
    elif 'splice_region_variant' in effects:
        return 'splice_region_variant'
    elif '5_prime_UTR_variant' in effects:
        return '5_prime_UTR_variant'
    elif '3_prime_UTR_variant' in effects:
        return '3_prime_UTR_variant'
    elif 'non_coding_transcript_exon_variant' in effects:
        return 'non_coding_transcript_exon_variant'
    elif 'intron_variant' in effects:
        return 'intron_variant'
    elif 'stop_retained_variant' in effects:
        return 'stop_retained_variant'
    elif 'synonymous_variant' in effects:
        return 'synonymous_variant'
    else:
        return '.'

# returns read depth, allele observation count, and percent allele observation count for all samples
def reads_stats(variants):
    stats = []
    
    for var in variants:
        DP = var.info['DP'] # read depth
        AOs = var.info['AO'] # allele observation count(s)

        # one allele associated with variant
        if len(AOs) == 1:
            AO = AOs[0]
            percent = round((AO/DP)*100, 2)
            stats.append([str(DP), str(AO), str(percent)])

        # more than one allele associated 
        elif len(AOs) > 1:
            percents = []
            
            # iterates through each allele 
            for AO in AOs:
                percent = round((AO/DP)*100, 2)
                percents.append(str(percent))
                
            AOs = [str(AO) for AO in AOs] 
            stats.append([str(DP), ','.join(AOs), ','.join(percents)])
    
    print('Total reads stats annotated.')        
    return stats

# returns allele frquency of each variant from ExAC API
def allele_freq(variants):
    AFs = []
    
    for var in variants:
        chrom = str(var.chrom)
        pos = str(var.pos)
        ref = var.ref
        alts = var.alts

        # one allele associated with variant
        if len(alts) == 1:

            # ExAC API request
            var_id = f'{chrom}-{pos}-{ref}-{alts[0]}'
            request = f'http://exac.hms.harvard.edu/rest/variant/variant/{var_id}'
            response = requests.get(request)
            text = json.loads(response.text)
            
            if 'allele_freq' in text.keys():
                AF = text['allele_freq']
                AFs.append(str(AF))   
            else:
                AFs.append('.')

        # more than one alleles associated with variant
        if len(alts) > 1:
            AFs_temp = []
            
            # iterate through each allele
            for alt in alts:

                # ExAC API request
                var_id = f'{chrom}-{pos}-{ref}-{alt}'
                request = f'http://exac.hms.harvard.edu/rest/variant/variant/{var_id}'
                response = requests.get(request)
                text = json.loads(response.text)
                
                if 'allele_freq' in text.keys():
                    AF = text['allele_freq']
                    AFs_temp.append(str(AF))   
                else:
                    AFs_temp.append('.')

            # concatenates elements of AFs_temp into one string and stores it in AFs
            AFs.append(','.join(AFs_temp))

    print('Allele frequencies annotated.')
    return AFs

# returns read depth, allele observation count, and percent allele observation count for each sample
def reads_per_sample_stats(variants):
    reads_per_sample_stats = []
    
    for var in variants:
        temp = []

        #iterates through each sample
        for sample in var.samples:
            DP = var.samples[sample]['DP'] # sample read depth
            AOs = var.samples[sample]['AO'] # sample allele observatopm count(s)

            # one allele associated with variant
            if len(AOs) == 1:
                AO = AOs[0]
                percent = round((AO/DP)*100, 2)
                temp.append([str(DP), str(AO), str(percent)])

            # more than one allele associated
            elif len(AOs) > 1:
                percents = []
                
                # iterates through each allele
                for AO in AOs:
                    percent = round((AO/DP)*100, 2)
                    percents.append(str(percent))
                    
                AOs = [str(AO) for AO in AOs] 
                temp.append([str(DP), ','.join(AOs), ','.join(percents)])

        # appends list of samples to reads_per_sample_stats
        reads_per_sample_stats.append(temp)

    print('Reads per sample stats annotated.')
    return reads_per_sample_stats

def main():

    # get parameters
    args = get_args()

    # load vcf data
    data = VariantFile(args.input)

    # convert vcf into python objects
    variants = get_variants(data)

    #get names of samples
    samples = get_samples(data)

    print(f'\nVCF file processed. There are {len(variants)} variants.\nBeginning annotation.\n')

    # run annotation functions
    var_types = variant_types(variants)
    consequences = effects(variants)
    reads = reads_stats(variants)
    AFs = allele_freq(variants)
    reads_per_sample = reads_per_sample_stats(variants)

    print('\nAnnotations completed. Writing tsv file.\n')

    # create header for tsv file
    columns = ['ref', 'alt', 'variant_type', 'effect','total_read_depth', 'total_alt_allele_reads', 
                'percent_total_alt_allele_reads(%)', 'allele_frequency']

    for sample in samples:
        columns.append(f'{sample}_read_depth')
        columns.append(f'{sample}_alt_allele_reads')
        columns.append(f'{sample}_percent_total_alt_allele_reads(%)')

    header = '\t'.join(columns)
    
    # writing results to output file
    output = open(args.output, 'w')
    output.write(header+'\n')

    index = range(len(variants))

    for i in index:
        row = '\t'.join([variants[i].ref, ','.join(variants[i].alts), var_types[i], consequences[i], \
              '\t'.join(reads[i]), AFs[i]]) + '\t'
    
        for j in range(len(samples)):
            row = row + '\t'.join(reads_per_sample[i][j]) + '\t'
            row.strip('\t')
    
        output.write(row + '\n')

    output.close()
    print('Finished.')


if __name__ == "__main__":
    main()