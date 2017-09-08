import docker
import os

def create_genesippr_container():
    """
    Starts the genesippr docker container
    """

    # Instantiate client to talk to Docker daemon
    client = docker.from_env()

    container = client.containers.run("genesipprv2:latest",
                                      detach=True,
                                      stdin_open=True,
                                      tty=True,
                                      volumes={'/mnt/nas': '/mnt/nas'})

    return container


def run_genesippr(container, input_directory, sequence_path, target_path):
    """
    Takes a genesippr container and runs the program
    """
    a = container.exec_run('python3 /geneSipprV2/sipprverse/genesippr/genesippr.py {0}' 
                           ' -s {1}' 
                           ' -t {2}' 
                           ' --detailedReports'.format(input_directory, sequence_path, target_path))
    print("{}".format(a.strip().decode('UTF-8')))

def mass_sipp(id_list):
    # Get list of folders to run genesippr on
    id_list = os.listdir('/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes')

    # Start the container
    container = create_genesippr_container()

    # Send commands to the container
    for id in id_list:
        print("\nRunning genesippr on {}...".format(id))
        run_genesippr(container=container,
                      input_directory='/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/{}'.format(id),
                      sequence_path='/mnt/nas/Forest/MG-RAST_Dataset_Analysis/metagenomes/{}'.format(id),
                      target_path='/mnt/nas/Forest/MG-RAST_Dataset_Analysis/db')