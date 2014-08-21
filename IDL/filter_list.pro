function filter_list,listfile
  spawn,'grep ">" '+listfile+" | grep -v # | awk '{print $2}' ",filters
  
end
