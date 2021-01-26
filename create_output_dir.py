def create_output_directory(filepath):

   '''
   Takes a filepath and finds the directory its in by stripping file off the end of the path.
   '''

   split = filepath.split("/")
   split.pop(-1)
   output_dir = ""
   for dir in split:
      output_dir+=dir+"/"
   return output_dir
