image_name=registry-vpc.cn-hangzhou.aliyuncs.com/rcecc/wgs
docker login --username=dggene --password=dg@12345 registry-vpc.cn-hangzhou.aliyuncs.com 
docker build -t $image_name:latest .
docker push $image_name:latest