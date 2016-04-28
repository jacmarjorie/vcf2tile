#!/usr/bin/env python

import subprocess, ConfigParser, json
import utils.helper as helper

class Loader:
  def __init__(self, config_file, debug=False):
    self.debug = debug
    self.readConfig(config_file)

  def readConfig(self, config_file):
    parser = ConfigParser.RawConfigParser()
    parser.read(config_file)
    self.NUM_PROCESSES    = parser.getint('loader', 'NUM_PROCESSES')
    self.EXEC             = parser.get('loader', 'EXECUTABLE')
    self.TILE_LOADER_JSON = parser.get('loader', 'TILE_LOADER_JSON')
    # If the load request was to a single tile db instance then
    # we can skip details on MPI, else MPI fields are required
    if self.NUM_PROCESSES == 1:
      return
    self.MPIRUN        = parser.get('mpi', 'MPIRUN')
    self.HOSTS         = parser.get('mpi', 'HOSTS')
    self.HOSTFLAG      = parser.get('mpi', 'HOSTFLAG')

    if parser.has_option('mpi', 'BTL_TCP_IF_INCLUDE'):
      self.IF_INCLUDE = parser.get('mpi', 'BTL_TCP_IF_INCLUDE')
    if parser.has_option('mpi', 'INCLUDE_ENV'):
      self.ENV      = parser.get('mpi', 'INCLUDE_ENV')

  def run(self, tile_loader_config_file):
    # if self.debug:
    #   helper.log("[Loader:Run] Starting mpirun subprocess")
    processArgs = list()
    if self.NUM_PROCESSES > 1:
      processArgs.extend([self.MPIRUN, "-np", str(self.NUM_PROCESSES),
                          self.HOSTFLAG, self.HOSTS])
      if self.IF_INCLUDE != None:
        processArgs.extend(["--mca", "btl_tcp_if_include", self.IF_INCLUDE])
      if self.ENV != None:
        for env in self.ENV.split(","):
          processArgs.extend(["-x", env])

    processArgs.extend([self.EXEC, tile_loader_config_file])
    # if self.debug:
    #   helper.log("Args: {0} ".format(processArgs))
    pipe = subprocess.Popen(processArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = pipe.communicate()

    if pipe.returncode != 0:
      raise Exception("subprocess run: {0}\nFailed with stdout: \n-- \n{1} \n--\nstderr: \n--\n{2} \n--".format(" ".join(processArgs), output, error))

def load2Tile(loader_config_file, callset_mapping_file, vid_mapping_file):
  loader = Loader(loader_config_file)

  with open(loader.TILE_LOADER_JSON, 'r') as loaderFP:
    loader_json_obj = json.load(loaderFP)
  loader_json_obj["callset_mapping_file"] = callset_mapping_file
  loader_json_obj["vid_mapping_file"]     = vid_mapping_file

  helper.writeJSON2File(loader_json_obj, loader.TILE_LOADER_JSON)
  print "Updated Tile Loader JSON file : {0}".format(loader.TILE_LOADER_JSON)
  loader.run(loader.TILE_LOADER_JSON)
  print "Tile DB loading complete"


def writeJSON2File(input_json, output_file):
  with open(output_file, "w") as outFP:
    json.dump(input_json, outFP, indent=2, separators=(',', ': '))

if __name__ == "__main__":
  import argparse 
  parser = argparse.ArgumentParser(description = "Load to tile db")

  parser.add_argument("-c", "--config", required=True, type=str,
                      help="input configuration file for invoking the tile loader")
  parser.add_argument("-d", "--debug", action="store_true", required=False,
                      help="debug mode flag")
  args = parser.parse_args()

  loader = Loader(args.config, args.debug)
  loader.run(loader.TILE_LOADER_JSON)
