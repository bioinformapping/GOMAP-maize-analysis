for f in $@
do
    echo $f
    out=`basename ${f%%.dot}`
    dot -Tpng -Gdpi=300 $f > $out.png
done
