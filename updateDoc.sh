git clone git@github.com:jmonlong/PopSV.git -b master ~/tmp
cp ~/tmp/*.md .
cp ~/tmp/public/* public/
rm -fr ~/tmp
mv README.md index.md

for f in *md
do
    cat $f | perl -ne '$_=~s/\]\(([0-9])/\]\({{ site\.baseurl }}$1/g; print $_;' > $f.tmp
    mv $f.tmp $f
done
