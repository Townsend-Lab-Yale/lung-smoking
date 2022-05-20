import synapseclient
import os
from locations import location_data

syn = synapseclient.Synapse()
username = input('Synapse username: ')
password = input('Synapse password: ')
syn.login(username, password)


# patient clinical data
entity = syn.get("syn24179660", downloadLocation= os.path.join(location_data, "genie_9"))
# sample data
entity = syn.get('syn24179661', downloadLocation= os.path.join(location_data, "genie_9"))
# MAF file
entity = syn.get('syn24179664', downloadLocation= os.path.join(location_data, "genie_9"))
# genomic information (genome coverage of TGS)
entity = syn.get('syn24179674', downloadLocation= os.path.join(location_data, "genie_9"))

