#update zoomx test as needed
tar -czvf zoomx.test.full.tgz -T zoomx.test.full.list
aws s3 cp zoomx.test.full.tgz s3://lixiabucket/zoomx.test.full.tgz

#update zoomx ami as needed
ssh -i "biolinux8_key_pair.pem" bkd@ec2-54-190-36-93.us-west-2.compute.amazonaws.com
#aws ec2 create-instance-export-task --instance-id i-0bcc981a6b1ab6744 --target-environment vmware --export-to-s3-task DiskImageFormat=vmdk,ContainerFormat=ova,S3Bucket=lixiabucket,S3Prefix=zoomxDemo.2017Apr26.
taskid=$(aws ec2 create-instance-export-task --instance-id i-0bcc981a6b1ab6744 --target-environment vmware --export-to-s3-task DiskImageFormat=vmdk,ContainerFormat=ova,S3Bucket=lixiabucket,S3Prefix=zoomxDemo. 2>&1 | grep -o export-i-[0-9a-z]*)
echo $taskid
aws ec2 describe-export-tasks --export-task-ids $taskid #show status
aws s3 cp s3://lixiabucket/zoomxDemo.$taskid.ova s3://lixiabucket/zoomxDemo.ova #working from MBP
