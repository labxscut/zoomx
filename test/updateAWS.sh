tar -czvf zoomx.test.full.tgz -T zoomx.test.full
aws s3 cp zoomx.test.full.tgz s3://lixiabucket/zoomx.test.full.tgz
aws ec2 create-instance-export-task --instance-id i-0bcc981a6b1ab6744 --target-environment vmware --export-to-s3-task DiskImageFormat=vmdk,ContainerFormat=ova,S3Bucket=lixiabucket,S3Prefix=zoomxDemo.2017Apr26.
#take a note of #"ExportTaskId": "export-i-fgciwtz3", "InstanceId": "i-0bcc981a6b1ab6744"
aws ec2 describe-export-tasks --export-task-ids export-i-fgciwtz3 #show status
aws s3 mv s3://lixiabucket/zoomxDemo.2017Apr26export-i-fgciwtz3.ova s3://lixiabucket/zoomxDemo.ova #working from MBP
