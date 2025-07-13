curl -o upload.sh https://raw.githubusercontent.com/juliansgottfried/Dengue/refs/heads/main/local_transfer/upload.sh
curl -o download.sh https://raw.githubusercontent.com/juliansgottfried/Dengue/refs/heads/main/local_transfer/download.sh
curl -o bind.R https://raw.githubusercontent.com/juliansgottfried/Dengue/refs/heads/main/local_transfer/helpers/bind.R

mkdir helpers fitting_folders
mv bind.R helpers

read -p "Enter your NYU HPC username: " user_name
echo "user_name=$user_name" > helpers/username.sh
