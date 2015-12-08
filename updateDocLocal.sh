cp ../../PopSV/*.md .
cp ../../PopSV/public/* public/
mv README.md index.md

for f in *md
do
    cat $f | perl -ne '$_=~s/\]\(([0-9])/\]\({{ site\.baseurl }}$1/g; print $_;' > $f.tmp
    mv $f.tmp $f
    cat $f | perl -ne '$_=~s/\]\(public/\]\({{ site\.baseurl }}public/g; print $_;' > $f.tmp
    mv $f.tmp $f
done
