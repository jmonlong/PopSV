cp ../../PopSV/*.md .
cp ../../PopSV/public/* public/
mv README.md index.md

for f in *md
do
    cat $f | perl -ne '$_=~s/\]\(([0-9])/\]\({{ site\.baseurl }}\/$1/g; print $_;' > $f.tmp
    cat $f.tmp | perl -ne '$_=~s/\]\(public/\]\({{ site\.baseurl }}\/public/g; print $_;' > $f
    rm $f.tmp
done

echo '---\nlayout: default\ntitle: PopSV\n---\n' > tmp.md
cat index.md >> tmp.md
mv tmp.md index.md

echo '---\nlayout: page\ntitle: Analysis steps\npermalink: /1-BasicWorkflow.md/\n---\n' > tmp.md
cat 1-BasicWorkflow.md >> tmp.md
mv tmp.md 1-BasicWorkflow.md

echo '---\nlayout: pagetoc\ntitle: Cluster management in R using BatchJobs\npermalink: /2-ClusterManagement-BatchJobs.md/\n---\n' > tmp.md
cat 2-ClusterManagement-BatchJobs.md >> tmp.md
mv tmp.md 2-ClusterManagement-BatchJobs.md

echo '---\nlayout: pagetoc\ntitle: Cluster management in R\npermalink: /2-ClusterManagement.md/\n---\n' > tmp.md
cat 2-ClusterManagement.md >> tmp.md
mv tmp.md 2-ClusterManagement.md

echo '---\nlayout: pagetoc\ntitle: Visualization\npermalink: /3-Visualization.md/\n---\n' > tmp.md
cat 3-Visualization.md >> tmp.md
mv tmp.md 3-Visualization.md

echo '---\nlayout: page\ntitle: Cancer analysis\npermalink: /4-Cancer.md/\n---\n' > tmp.md
cat 4-Cancer.md >> tmp.md
mv tmp.md 4-Cancer.md

echo '---\nlayout: pagetoc\ntitle: Frequently Asked Questions\npermalink: /5-FAQ.md/\n---\n' > tmp.md
cat 5-FAQ.md >> tmp.md
mv tmp.md 5-FAQ.md

echo '---\nlayout: page\ntitle: Publications\npermalink: /6-Publication.md/\n---\n' > tmp.md
cat 6-Publication.md >> tmp.md
mv tmp.md 6-Publication.md

echo '---\nlayout: page\ntitle: About\npermalink: /9-About.md/\n---\n' > tmp.md
cat 9-About.md >> tmp.md
mv tmp.md 9-About.md

