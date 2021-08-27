import synapseclient
syn = synapseclient.Synapse()
username = input('Synapse username: ')
password = input('Synapse password: ')
syn.login(username, password)
entity = syn.get("syn24179660", downloadLocation="../data/genie_9")
entity = syn.get('syn24179661', downloadLocation="../data/genie_9")
entity = syn.get('syn24179664', downloadLocation="../data/genie_9")
entity = syn.get('syn24179674', downloadLocation="../data/genie_9")
