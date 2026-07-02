dir="./data/illumina/05757/2010_test"
out="$(basename "$(dirname "$dir")")_$(basename "$dir").tsv"
ids=$(cd "$dir" && ls *_1.fastq | sed 's/_1\.fastq$//' | sort)
prefix=$(printf "%s\n" "$ids" | awk '
NR==1 { p=$0; next }
{
    while (index($0,p)!=1) p=substr(p,1,length(p)-1)
}
END { print p }
')
prefix="${prefix%_}"
echo "[INFO] Common prefix: $prefix"


if [ -z "$prefix" ]; then
    echo "[ERROR] No common prefix found between FASTQ names." >&2
    exit 1
fi

(cd "$dir" && printf "%s\n" *_1.fastq) | sort | awk -v sample="$prefix" '
BEGIN { OFS="\t"; fc=1 }
{
    r1=$0
    lib=r1; sub(/_1\.fastq$/, "", lib)
    r2=lib"_2.fastq"
    print sample, lib, r1, r2, "ILLUMINA", sprintf("FC%03d",fc), 1, sprintf("BC%02d",fc), "NA", "RUN"fc, "SEQ"fc
    fc++
}' > "$out"