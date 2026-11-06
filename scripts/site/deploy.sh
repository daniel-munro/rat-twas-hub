set -e

rsync -av --delete --exclude='data' jekyll/_site/ root@ratgtex.org:/var/www/twas.ratgtex.org/ $1
rsync -av --delete jekyll/data/ root@ratgtex.org:/mnt/blockstorage/twas/data/ $1
ssh root@ratgtex.org "ln -sf /mnt/blockstorage/twas/data /var/www/twas.ratgtex.org/data"
ssh root@ratgtex.org "systemctl restart nginx"
