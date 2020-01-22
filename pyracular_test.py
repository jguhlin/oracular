import pyracular
NtFasta = pyracular.NtFasta(17, "/mnt/data/nt/nt.gz", 5, 1024, 8192, 8192, 2)
for i in range(100000):
	batch = NtFasta.next_batch();
	if len(batch) == 0:
		print("Finished!")
		break
	print(batch[0]["ID"])

