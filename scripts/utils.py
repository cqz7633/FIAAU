import os 
import time

def check_path(path):
    if os.path.isabs(path):
        path = path
    else:
        path = os.getcwd()  + "/" + path
    return path

def create_file(input_path):
    if not os.path.exists(input_path):
        os.mkdir(input_path)

def bam2bw(bam, bw, p, bin=1):
    commond = "bamCoverage -p %d --bam %s -o %s --binSize %d" % (p, bam, bw, bin)
    return commond

def bw2wig(bw,wig):
    commond = "bigWigToWig %s %s" % (bw, wig)
    return commond

def table_cond(table_path):
    bam_table_f = open(table_path, "r")
    bam_table_list = [i.strip() for i in bam_table_f.readlines() if i.strip() != ""]
    bam_table_f.close()
    cond_dict = {}
    for i in bam_table_list:
        tmp = i.split("\t")
        bam = tmp[0]
        type = tmp[1]
        if type not in cond_dict:
            cond_dict[type] = []
        cond_dict[type].append(bam)
    return  cond_dict

def check_split(key, paired):
    if paired == "y":
        if "_R1" in key:
            split_key = "_R1"
        elif ".R1" in key:
            split_key = ".R1"
        # elif "_1" in key:
        #     split_key = "_1"
        elif ".fq" in key:
            split_key = ".fq"
        elif ".fast" in key:
            split_key = ".fast"
        elif ".fa" in key:
            split_key = ".fa"
        return split_key
    if paired == "n":
        if ".fq" in key:
            split_key = ".fq"
        elif ".fast" in key:
            split_key = ".fast"
        elif ".fa" in key:
            split_key = ".fa"
        return split_key
    
def write_command(file, command):
    file.write("echo " + '"' + command + '"' + "\n")
    file.write(command + "\n")

def bam2bdg(geno_size, cond_list, bam_dict, out_path):
    bdg_dict = {}
    com_list = []
    for cond in cond_list:
        bam_list = bam_dict[cond]
        if cond not in bdg_dict:
            bdg_dict[cond] = []
        for bam in bam_list:
            tmp = bam.split("/")[-1]
            prefix = tmp.split(".bam")[0]
            bdg = out_path + prefix + ".bdg"
            bdg_dict[cond].append(bdg)
            com = "genomeCoverageBed -bg -ibam %s -g %s -split > %s" % (bam, geno_size, bdg)
            com_list.append(com)
    return bdg_dict, com_list


def write_dapars_config(config_path, bed, cond1, cond2, out_path, out_prefix, cov_cutoff):
    config_list = ["Annotated_3UTR=%s\n" % bed, \
            "Group1_Tophat_aligned_Wig=%s\n" % cond2, \
            "Group2_Tophat_aligned_Wig=%s\n" % cond1, \
            "Output_directory=%s\n" % out_path, \
            "Output_result_file=%s\n" % out_prefix, \
            "\n", \
            "#Parameters\n", \
            "Num_least_in_group1=1\n", \
            "Num_least_in_group2=1\n", \
            "Coverage_cutoff=%d\n" % cov_cutoff, \
            "FDR_cutoff=0.05\n", \
            "PDUI_cutoff=0\n", \
            "Fold_change_cutoff=0\n"]
    config_path.writelines(config_list)

def dapars_process(dapars_main, cond_list, bam_dict, out_path, core, dapars_bed, cov_cutoff, bin_size):
    dapars_file = out_path + "/DaPars/"
    dapars_data = dapars_file + "process/"
    dapars_data_bam = dapars_data + "bam/"
    dapars_data_bw = dapars_data + "bw/"
    dapars_data_wig = dapars_data + "/wig/"
    dapars_out = dapars_file + "out/"
    dapars_sh = dapars_file + "dapars.sh"
    dapars_config = dapars_file + "dapars_config.txt"
    
    create_file(dapars_file)
    create_file(dapars_data)
    create_file(dapars_data_bam)
    create_file(dapars_data_bw)
    create_file(dapars_data_wig)
    create_file(dapars_out)
    dapars_f = open(dapars_sh, "w")
    write_command(dapars_f,"cd %s" % dapars_file)
    
    bam_dict2 = {}
    ### ln bam to out dir
    for i in cond_list:
        if i not in bam_dict2:
            bam_dict2[i] = []
        for j in bam_dict[i]:
            tmp = j.split("/")[-1]
            ln_out = "%s%s" % (dapars_data_bam, tmp)
            ln_com = "ln -s %s %s" % (j, ln_out)
            write_command(dapars_f,ln_com)
            bam_dict2[i].append(ln_out)
    
    num = 1
    cond1_list = []
    cond2_list = []
    for i in cond_list:
        for j in bam_dict2[i]:
            bai_path = j + ".bai"
            if not os.path.exists(bai_path):
                bamindex_com = "samtools index -@ %d %s %s" % (core, j, bai_path)
                write_command(dapars_f,bamindex_com)
            tmp = j.split("/")[-1]
            prefix = tmp.split(".bam")[0]
            bw = dapars_data_bw + prefix + ".bw"
            wig = dapars_data_wig + prefix + ".wig"
            if not os.path.exists(bw):
                bam2bw_com = bam2bw(j, bw, core, bin_size)
                write_command(dapars_f,bam2bw_com)
            if not os.path.exists(wig):
                bw2wig_com = bw2wig(bw,wig)
                write_command(dapars_f,bw2wig_com)
            if num == 1:
                cond1_list.append(wig)
            elif num == 2:
                cond2_list.append(wig)
        if num == 1:
            cond1 = ",".join(cond1_list)
        elif num == 2:
            cond2 = ",".join(cond2_list)
        num += 1
    out_prefix = "%s_vs_%s" % (cond_list[1], cond_list[0])
    print(out_prefix)
    dapars_config_f = open(dapars_config, "w")
    write_dapars_config(dapars_config_f, dapars_bed, cond1, cond2, dapars_out, out_prefix, cov_cutoff)
    dapars_config_f.close()
    dapars_run = "python  %s %s" % (dapars_main, dapars_config)
    write_command(dapars_f, dapars_run)
    sub_dapars_com = "nohup sh %s > %sdapars.log 2>&1 &" % (dapars_sh, dapars_file)
    dapars_f.close()
    dapars_res = dapars_out + out_prefix + "_All_Prediction_Results.txt"
    return sub_dapars_com, dapars_res

def csiutr_process(cond_list, bam_dict, out_path, csiutr_arg, csiutr_bed, csiutr_anno, csiutr_anno_dir):
    csiutr_file = out_path + "/CSI_UTR/"
    csiutr_data = csiutr_file + "process/"
    csiutr_out = csiutr_file + "out/"
    create_file(csiutr_file)
    create_file(csiutr_data)
    create_file(csiutr_out)
    csiutr_sh = csiutr_file + "csiutr.sh"
    csi_bam = csiutr_file + "csi_bam.txt"
        
    csiutr_sh_f = open(csiutr_sh, "w")
    write_command(csiutr_sh_f,"cd %s" % csiutr_file)
    csi_bam_f = open(csi_bam, "w")
    for i in cond_list:
        rep_num = 1
        for j in bam_dict[i]:
            tmp = j.split("/")[-1]
            prefix = tmp.split(".bam")[0]
            ln_out = "%s%s" % (csiutr_data, tmp)
            ln_com = "ln -s %s %s" % (j, ln_out)
            bam_line = "%s\t%s\t%d" % (prefix, i, rep_num)
            write_command(csiutr_sh_f,ln_com)
            csi_bam_f.write(bam_line + "\n")
            rep_num += 1
    ln_anno_com = "ln -sf %s %s/annotations/Mm10.CSIs.annot.bed" % (csiutr_anno, csiutr_anno_dir)
    write_command(csiutr_sh_f, ln_anno_com)
    ln_bed_com = "ln -sf %s %s/locations/Mm10.CSIs.bed" % (csiutr_bed, csiutr_anno_dir)
    write_command(csiutr_sh_f, ln_bed_com)
    csi_run = "CSI-UTR %s -sample_info %s -bed %s -annot %s -out %s -data_dir %s" % (csiutr_arg, csi_bam, csiutr_bed, csiutr_anno, csiutr_out, csiutr_data)
    write_command(csiutr_sh_f, csi_run)
    csiutr_sh_f.close()
    sub_csiutr_com = "nohup sh %s > %scsiutr.log 2>&1 &" % (csiutr_sh, csiutr_file)
    csiutr_res = csiutr_out + "DifferentialExpression/WITHIN_UTR/" + cond_list[1] + "_VS_" + cond_list[0] + "_CSI.WITHINUTR.diff.txt"
    return sub_csiutr_com, csiutr_res

def diffutr_process(script, cond_list, bam_table, out_path, anno, method, cov_cutoff, reads_len):
    method_dir = out_path  + "/"+ method +"/"
    method_data = method_dir + "process/"
    method_out = method_dir + "out/"
    create_file(method_dir)
    create_file(method_data)
    create_file(method_out)
    method_sh = method_dir + method + ".sh"
    method_sh_f = open(method_sh, "w")
    write_command(method_sh_f,"cd %s" % method_dir)
    run_script_com = "Rscript %s -b %s -o %s -a %s -c %d -r %s" % (script, bam_table, method_dir, anno, cov_cutoff, reads_len)
    write_command(method_sh_f,run_script_com)
    method_sh_f.close()
    sub_script_com = "nohup sh %s > %s%s.log 2>&1 &" % (method_sh, method_dir, method)
    diffutr_res = method_out + cond_list[1] + "_vs_" + cond_list[0] + ".txt"
    return sub_script_com, diffutr_res

def qapa_process(fq_table, cond_list, out_path, paired, qapa_genome_fa, qapa_utr, core, qapa_ident, qapa_diff_script):
    qapa_file = out_path + "/QAPA/"
    qapa_data = qapa_file + "process/"
    qapa_out = qapa_file + "out/"
    create_file(qapa_file)
    create_file(qapa_data)
    create_file(qapa_out)
    qapa_sh = qapa_file + "qapa.sh"
    prefix_table = qapa_file + "prefix_table.txt"
    
    qapa_sh_f = open(qapa_sh, "w")
    prefix_table_f = open(prefix_table, "w")
    fq_table_f = open(fq_table, "r")
    qapa_sh_f.write("#!/bin/bash\n")
    write_command(qapa_sh_f,"cd %s" % qapa_file)
    ###extract 3utr seq
    qapa_utr_fa = "%s/utr_sequences.fa" % qapa_data
    utr_fa_com = "qapa fasta -f %s %s %s" % (qapa_genome_fa, qapa_utr, qapa_utr_fa)
    write_command(qapa_sh_f,utr_fa_com)
    
    ###crate 3utr salmon index
    qapa_index = "%s/3utr_salmon" % (qapa_data)
    index_com = "salmon index -t %s -i %s -p %d" % (qapa_utr_fa, qapa_index, core)
    write_command(qapa_sh_f,index_com)
    
    cond_prefix_dict = {}
    fq_table_f = open(fq_table, "r")
    fq_table_list = [i.strip() for i in fq_table_f.readlines() if i.strip() != ""]
    fq_table_f.close()
    n = 1
    for line1 in fq_table_list:
        if paired != "y":
            tmp = line1.split("\t")
            fq = tmp[0]
            fq_file = fq.split("/")[-1] #R1 file name split by full path
            split_key = check_split(fq_file, paired) #find split key of prefix
            prefix = fq_file.split(split_key)[0] #prefix by R1
            cond = tmp[1]
            prefix_table_f.write("%s\t%s\n" % (prefix, cond))
            if cond not in cond_prefix_dict:
                cond_prefix_dict[cond] = []
            cond_prefix_dict[cond].append(prefix)
            #rna-seq command
            salmon = "salmon quant -i %s -l A -r %s --validateMappings  -o %s%s -p %d" % (qapa_index, fq, qapa_data, prefix, core)
            #write rna-seq to the sh file
            if not os.path.exists(qapa_data+prefix+"/quant.sf"):
                write_command(qapa_sh_f, salmon)
                qapa_sh_f.write("echo '________%s processed end________'\n\n" % prefix)
        
        #pair end processed
        elif paired == "y":
            tmp = line1.split("\t")
            #R1 file name
            if n % 2 == 1:
                fq1 = tmp[0]
            #R1 and R2 are all find
            elif n % 2 == 0:
                fq2 = tmp[0]
                prefix = fq1.split("/")[-1] #R1 file name split by full path
                split_key = check_split(prefix, paired) #find split key of prefix
                prefix = prefix.split(split_key)[0] #prefix by R1
                cond = tmp[1]
                prefix_table_f.write("%s\t%s\n" % (prefix, cond))
                if cond not in cond_prefix_dict:
                    cond_prefix_dict[cond] = []
                cond_prefix_dict[cond].append(prefix)
                #rna-seq command
                salmon = "salmon quant -i %s -l A -1 %s -2 %s --validateMappings  -o %s%s -p %d" % (qapa_index, fq1, fq2, qapa_data, prefix, core)
                if not os.path.exists(qapa_data+prefix+"/quant.sf"):
                    write_command(qapa_sh_f, salmon)
                    qapa_sh_f.write("echo '________%s processed end________'\n\n" % prefix)
        n += 1
    prefix_table_f.close()
    qapa_compare = "%s_vs_%s" % (cond_list[1], cond_list[0])
    qapa_pau = "%s_pau.txt" % (qapa_compare)
    run_qapa_com = "qapa quant --db %s %s*/quant.sf > %s%s" % (qapa_ident, qapa_data, qapa_out, qapa_pau)
    run_qapa_test_com = "Rscript %s -i %s%s -p %s -c %s -o %s" % (qapa_diff_script, qapa_out, qapa_pau, prefix_table, qapa_compare, qapa_out)
    write_command(qapa_sh_f, run_qapa_com)
    write_command(qapa_sh_f, run_qapa_test_com)
    qapa_sh_f.close()
    sub_qapa_sh = "nohup bash %s > %sqapa.log 2>&1 &" % (qapa_sh, qapa_file)
    qapa_res = qapa_out + qapa_compare + "_pau_result.txt"
    return sub_qapa_sh, qapa_res

def apatrap_process(cond_list, bam_dict, out_path, identify_script, predict_script, deapa_script, geno_size, apatrap_genemodel, cov_cutoff):
    apatrap_file = out_path + "/APAtrap/"
    apatrap_data = apatrap_file + "process/"
    apatrap_out = apatrap_file + "out/"
    create_file(apatrap_file)
    create_file(apatrap_data)
    create_file(apatrap_out)
    apatrap_sh = apatrap_file + "apatrap.sh"
    apatrap_sh_f = open(apatrap_sh, "w")
    write_command(apatrap_sh_f, "cd %s" % apatrap_file)
    cov_cutoff = int(cov_cutoff * 10)
    if cov_cutoff < 10:
        cov_cutoff = 10
    ### bam to bedgraph 
    bdg_dict, com_list = bam2bdg(geno_size, cond_list, bam_dict, apatrap_data)
    write_command(apatrap_sh_f, "\n".join(com_list))

    # 1.identifyDistal3UTR
    bdg_input = " ".join(bdg_dict[cond_list[0]]+bdg_dict[cond_list[1]])
    identy_out = apatrap_out + "3utr.refine.bed"
    identy_com = "%s -i %s -m %s -o %s -c 0.05" % (identify_script, bdg_input, apatrap_genemodel, identy_out)
    write_command(apatrap_sh_f, identy_com)
    
    # 2.predict APA
    cond1_num = len(bdg_dict[cond_list[0]])
    cond2_num = len(bdg_dict[cond_list[1]])
    bdg_cond1_list = bdg_dict[cond_list[0]]
    bdg_cond2_list = bdg_dict[cond_list[1]]
    bdg_cond1_input = " ".join(bdg_cond1_list)
    bdg_cond2_input = " ".join(bdg_cond2_list)
    predictapa_out = apatrap_out + "%s_vs_%s_predictAPA.txt" % (cond_list[1],cond_list[0])
    ###use refind utr
    predictapa_com = "%s  -i %s %s -g 2 -n %d %d -u %s -o %s -c %d" % (predict_script, bdg_cond1_input, bdg_cond2_input, cond2_num, cond1_num, identy_out, predictapa_out, cov_cutoff)
    write_command(apatrap_sh_f, predictapa_com)
    
    # 3. deAPA
    deapa_out = apatrap_out + "%s_vs_%s_deAPA.txt" % (cond_list[1],cond_list[0])
    deapa_com = "Rscript %s -i %s -o %s -c %d" % ( deapa_script, predictapa_out, deapa_out, cov_cutoff)
    write_command(apatrap_sh_f, deapa_com)
    
    sub_apatrap_sh = "nohup bash %s > %sapatrap.log 2>&1 &" % (apatrap_sh, apatrap_file)

    return sub_apatrap_sh, deapa_out

def labrat_process(fq_table, cond_list, out_path, labrat_fa, labrat_gff, labrat_db, labrat_seq, paired, maketf, core):
    labrat_file = out_path + "/LABRAT/"
    labrat_data = labrat_file + "process/"
    labrat_out = labrat_file + "out/"
    create_file(labrat_file)
    create_file(labrat_data)
    create_file(labrat_out)
    labrat_sh = labrat_file + "labrat.sh"
    sample_info = labrat_file + "sample_info.txt"
    labrat_sh_f = open(labrat_sh,"w")
    write_command(labrat_sh_f,"cd %s" % labrat_data)
    labrat_new_gff = "%s/%s" % (labrat_data, os.path.basename(labrat_gff))
    ls_gff_com = "ln -s %s %s" % (labrat_gff, labrat_new_gff)
    write_command(labrat_sh_f, ls_gff_com)

    ###fq dict
    fq1_list = []
    fq2_list = []
    sample_name_list = []
    flag_list = []
    fq_table_f = open(fq_table, "r")
    sample_info_f = open(sample_info,"w")
    sample_info_f.write("sample\tcondition\t\n")
    n = 1
    for line in fq_table_f:
        tmp = line.strip().split("\t")
        fq = tmp[0]
        cond = tmp[1]
        fq_name = fq.split("/")[-1]
        ln_out = "%s/%s" % (labrat_data, fq_name)
        ln_com = "ln -s %s %s" % (fq, ln_out)
        write_command(labrat_sh_f, ln_com)
        if paired == "y":
            if n % 2 != 0:
                fq1_list.append(fq_name)
                if cond not in flag_list:
                    sample_num = 1
                    flag_list.append(cond)
                sample = "%s_%d" % (cond, sample_num)
                sample_name_list.append(sample)
                sample_out = "%s\t%s" % (sample,cond)
                sample_info_f.write(sample_out + "\n")
                sample_num += 1
            else:
                fq2_list.append(fq_name)
        else:
            fq1_list.append(fq_name)
            if cond not in flag_list:
                sample_num = 1
                flag_list.append(cond)
            sample = "%s_%d" % (cond, sample_num)
            sample_name_list.append(sample)
            sample_out = "%s\t%s" % (sample,cond)
            sample_info_f.write(sample_out + "\n")
            sample_num += 1
        n += 1
    fq_table_f.close()
    sample_info_f.close()
    read1_fq = ",".join(fq1_list)
    read2_fq = ",".join(fq2_list)
    sample_name_combine = ",".join(sample_name_list)
    mkfa_com = "LABRAT.py --mode makeTFfasta --gff %s --genomefasta %s --lasttwoexons --librarytype RNAseq --threads %d" \
               % (labrat_new_gff, labrat_fa, core)
    txfa = "%s/TFseqs.fasta" % (labrat_data)
    if maketf == "y":
        if paired == "y":
            runsalm_com = "LABRAT.py --mode runSalmon --librarytype RNAseq --txfasta %s --reads1 %s --reads2 %s --samplename %s --threads %d" \
                        % (txfa, read1_fq, read2_fq, sample_name_combine, core)
        else:
            runsalm_com = "LABRAT.py --mode runSalmon --librarytype RNAseq --txfasta %s --reads1 %s --samplename %s --threads %d" \
                        % (txfa, read1_fq, sample_name_combine, core)
    elif maketf == "n":
        cp_db_com = "cp %s %s" % (labrat_db, labrat_data)
        if paired == "y":
            runsalm_com = "LABRAT.py --mode runSalmon --librarytype RNAseq --txfasta %s --reads1 %s --reads2 %s --samplename %s --threads %d" \
                        % (labrat_seq, read1_fq, read2_fq, sample_name_combine, core)
        else:
            runsalm_com = "LABRAT.py --mode runSalmon --librarytype RNAseq --txfasta %s --reads1 %s --samplename %s --threads %d" \
                        % (labrat_seq, read1_fq, sample_name_combine, core)
    salmon_dir = labrat_data
    calpsi_com = "LABRAT.py --mode calculatepsi --gff %s --librarytype RNAseq  --salmondir %s --sampconds %s --conditionA %s --conditionB %s" \
                % (labrat_new_gff, salmon_dir, sample_info, cond_list[0], cond_list[1])
    labrat_res = "%s/%s_vs_%s.LABRAT.psis.pval" % (labrat_out, cond_list[1], cond_list[0])
    move_res_com = "mv %s/LABRAT.psis.pval %s" % (labrat_data, labrat_res)
    if maketf == "y":
        write_command(labrat_sh_f, mkfa_com)
    elif maketf == "n":
        write_command(labrat_sh_f, cp_db_com)
    write_command(labrat_sh_f, runsalm_com)
    write_command(labrat_sh_f, calpsi_com)
    write_command(labrat_sh_f, move_res_com)
    labrat_sh_f.close()
    sub_labrat_sh = "nohup bash %s > %slabrat.log 2>&1 &" % (labrat_sh, labrat_file)
    return sub_labrat_sh, labrat_res
