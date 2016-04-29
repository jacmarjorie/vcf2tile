import sys, os, vcf, uuid, json, subprocess

from collections import OrderedDict

CONST_TILEDB_FIELDS = OrderedDict()
CONST_TILEDB_FIELDS["END"]             = { "vcf_field_class" : ["INFO"],          "type": "int" }
CONST_TILEDB_FIELDS["BaseQRankSum"]    = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["ClippingRankSum"] = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["MQRankSum"]       = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["ReadPosRankSum"]  = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["MQ"]              = { "vcf_field_class" : ["INFO"],          "type":"float" }
CONST_TILEDB_FIELDS["MQ0"]             = { "vcf_field_class" : ["INFO"],          "type":"int"   }
CONST_TILEDB_FIELDS["AF"]              = { "vcf_field_class" : ["INFO"],          "type":"float", "length":"A" }
CONST_TILEDB_FIELDS["AN"]              = { "vcf_field_class" : ["INFO"],          "type":"int",   "length":1 }
CONST_TILEDB_FIELDS["AC"]              = { "vcf_field_class" : ["INFO"],          "type":"int",   "length":"A" }
CONST_TILEDB_FIELDS["DP"]              = { "vcf_field_class" : ["INFO","FORMAT"], "type":"int"   }
CONST_TILEDB_FIELDS["MIN_DP"]          = { "vcf_field_class" : ["FORMAT"],        "type":"int"   }
CONST_TILEDB_FIELDS["GQ"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int"   }
CONST_TILEDB_FIELDS["SB"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":4 }
CONST_TILEDB_FIELDS["AD"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":"R" }
CONST_TILEDB_FIELDS["PL"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":"G" }
CONST_TILEDB_FIELDS["GT"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":"P" }
CONST_TILEDB_FIELDS["PS"]              = { "vcf_field_class" : ["FORMAT"],        "type":"int",   "length":1 }


class VCF:

  def __init__(self, filename, loader_config):
    """
    init uses reads in the vcf
    """
    self.filename = os.path.abspath(filename)

    # config information
    with open(loader_config) as conf:
      config = json.load(conf)
    
    self.callset_map = config['callset_mapping_file']
    self.vid_map = config.get('vid_mapping_file', None)
    # assuming only one tiledb instance
    self.workspace = config['column_partitions'][0]['workspace']
    self.array = config['column_partitions'][0]['array'] 

    # TileDB loader specific 
    self.callset_mapping = dict()

    # just read callsets from header by default
    self.sampleTag = False


  def __enter__(self):

    self.file = open(self.filename, 'rb')
    self.reader = vcf.Reader(self.file)

    return self

  def __exit__(self, exc_type, exc_value, traceback):
    self.file.close()


  def writeVIDMappingFile(self, fields_dict=CONST_TILEDB_FIELDS):
      vid_mapping = OrderedDict()
      vid_mapping["fields"] = fields_dict
      vid_mapping["contigs"] = OrderedDict()
      contigs = vid_mapping["contigs"]
      offset = 0
      # use vcf header to construct contigs
      for reference in self.reader.contigs:
        contigs[reference] = {"length": self.reader.contigs[reference].length,
                                   "tiledb_column_offset": offset}
        offset += (long(self.reader.contigs[reference].length) + 1000)

      writeJSON2File(vid_mapping, self.vid_map)


  def getCallSets(self, callset_check, row_counter=0):

    callsets = OrderedDict()

    for i in range(0, len(self.reader.samples)):
      if self.sampleTag is False:
        name = self.reader.samples[i]
        idx_in_file = i
      else:
        name = self.reader.metadata['SAMPLE'][i]['SampleName']
        idx_in_file = self.reader.samples.index(self.reader.metadata['SAMPLE'][i]['ID'])

      # check duplicate
      if name in callset_check:
        nuuid = str(uuid.uuid4())
        sys.stderr.write('Duplicate callset name '+name+' : appending _'+nuuid+'\n');
        name += ('_'+nuuid)

      callsets[name] = OrderedDict()
      callsets[name]['row_idx'] = row_counter
      callsets[name]['idx_in_file'] = idx_in_file
      callsets[name]['filename'] = self.filename 
      print callsets[name]
      print 'callset', name, 'with', str(row_counter)

      row_counter += 1


    return callsets, row_counter


def writeJSON2File(input_json, output_file):
  with open(output_file, "w") as outFP:
    json.dump(input_json, outFP, indent=2, separators=(',', ': '))


if __name__ == "__main__":
  import argparse 
  parser = argparse.ArgumentParser(description = "Load a list of VCFs into TileDB.")

  # includes path to vid_map and callset_map
  # workspace, array
  parser.add_argument("-l", "--loader", required=True, type=str,
                      help="base loader file")
  parser.add_argument("-i", "--inputs", nargs='+', type=str, required=True, help="VCF files to be imported.")
  parser.add_argument("-s", "--sampleTag", action="store_true", required=False, help="Use SAMPLE tag in VCF header to name callsets.")
  
  args = parser.parse_args()

  callset_mapping = {}
  callset_mapping['callsets'] = OrderedDict()
  callsets = callset_mapping['callsets']

  counter = 0
  rc = 0
  for input_file in args.inputs:

    with VCF(input_file, args.loader) as vc:
      vc.sampleTag = args.sampleTag
      # write vid mapping file for the first
      # ideally one array should have one vid array associated with it
      if counter == 0:
        vc.writeVIDMappingFile()
        callset_map_file = vc.callset_map

      new_callset, rc = vc.getCallSets(callsets, row_counter=rc)
      callsets.update(new_callset)
      counter += 1
      
  writeJSON2File(callset_mapping, callset_map_file)

  processArgs = ['vcf2tiledb', os.path.abspath(args.loader)]
  # load the files
  pipe = subprocess.Popen(processArgs, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  output, error = pipe.communicate()

  if pipe.returncode != 0:
    raise Exception("subprocess run: {0}\nFailed with stdout: \n-- \n{1} \n--\nstderr: \n--\n{2} \n--".format(" ".join(processArgs), output, error))


