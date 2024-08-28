git clone https://github.com/Illumina/SpliceAI.git
cd SpliceAI
python setup.py install


spliceai -I /disk1/wenqing/SpliceAI/examples/input.vcf -O ./output.vcf -R /disk1/wenqing/tmp_data/hg19/hg19.fa -A /disk1/wenqing/SpliceAI/spliceai/annotations/grch37.txt


####缺少tenserrt库
tar -xzvf TensorRT-8.6.1.6.Linux.x86_64-gnu.cuda-11.8.tar.gz # 解压文件 
# 将lib添加到环境变量里面 
vim ~/.bashrc 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./TensorRT-8.6.1.6/lib 
source ~/.bashrc 
 
https://developer.download.nvidia.com/compute/cuda/10.2/Prod/local_installers/cuda_10.2.89_440.33.01_linux.run
sudo sh cuda_10.2.89_440.33.01_linux.run


# 或 直接将 TensorRT-8.6.1.6/lib 添加到 cuda/lib64 里面 
cp -r ./lib/* /usr/local/cuda/lib64/ 
 
# 安装python的包 
cd TensorRT-8.6.1.6/python 
pip install tensorrt-xxx-none-linux_x86_64.whl




curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh