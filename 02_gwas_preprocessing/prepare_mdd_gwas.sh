# paths
RAW="path/to/pgc-mdd2025_no23andMe_eur.tsv.gz"
OUT="path/to/PGC_MDD2025_no23andMe_EUR_ready.tsv"

# prepare GWAS
gunzip -c "${RAW}" \
| awk 'BEGIN{FS=OFS="\t"}
  /^##/ {next}
  /^#CHROM/ {
    gsub(/^#/, "", $1)
    for(i=1;i<=NF;i++) pos[$i]=i
    req[1]="ID"; req[2]="EA"; req[3]="NEA"; req[4]="BETA"; req[5]="SE"; req[6]="PVAL"; req[7]="NEFF"
    for(j=1;j<=7;j++){
      if(!(req[j] in pos)){
        print "ERR: missing column " req[j] > "/dev/stderr"
        exit 1
      }
    }
    print "SNP","A1","A2","BETA","SE","P","N","Z"
    next
  }
  {
    snp=$(pos["ID"]); a1=$(pos["EA"]); a2=$(pos["NEA"])
    beta=$(pos["BETA"]); se=$(pos["SE"]); p=$(pos["PVAL"]); n=$(pos["NEFF"])

    if(snp=="" || snp==".") next
    if(a1=="" || a2=="") next
    if(beta=="" || se=="" || se=="0") next
    if(p=="" || p=="NA") next
    if(n=="" || n=="NA") next

    z = beta / se
    print snp, a1, a2, beta, se, p, n, z
  }' > "${OUT}"