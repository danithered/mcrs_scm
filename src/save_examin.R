# Load the library xml2
library(xml2)# Read the xml file

#Functions
get_child <- function(node, name, as = NA){
  w <- which(name == xml_name(xml_children(node)))
  if (length(w) != 1) warning("get_child: number of correct children does ot equals to 1!")
  if(!is.na(as)) {
    if(as == "int") return(xml_integer(xml_child(node, w)))
    else if(as == "double") return(xml_double(xml_child(node, w)))
    else if(as == "string") return(xml_text(xml_child(node, w)))
    else if(as == "logical") return( ifelse(xml_integer(xml_child(node, w)) == 0, F, T ))
  }
  return(xml_child(node, w))
}

children <- function(node){
  xml_name(xml_children(node))
}

# Get data

data= read_xml('/home/danielred/data/programs/mcrs_to_scm/OUT/test/SAVE/c0.xml')
d <- xml_child(data, 1) # mooving to mcrscm

time <- get_child(d, "time", "int")
size <- get_child(d, "sim.size", "int")

cells <- get_child(d, "cells")

# get table
table <- do.call(rbind, lapply(1:xml_length(cells), function(no_cell){
  cell <- xml_child(cells, no_cell)

  list( alive = get_child(cell, "cell.alive", "logical"),
        leftover = get_child(cell, "cell.leftover", "double"),
        M = get_child(cell, "metabolism", "double"),
        no_reps = get_child(cell, "cell.reps") |> get_child("count", "int"),
        reps = get_child(cell, "cell.reps")
  )
}))

# examine reps
reps = table[1,]$reps

rep_table = do.call(rbind, lapply(which(children(reps) == "item"), function(rnum){
  rep <- xml_child(reps, rnum)
  out <- list(seq = get_child(rep, "repl.seq", "string"),
       str = get_child(rep, "str", "string"),
       mfe = get_child(rep, "repl.mfe", "double"),
       Pfold = get_child(rep, "repl.Pfold", "double"),
       Pdeg = get_child(rep, "repl.Pdeg", "double"),
       R = get_child(rep, "repl.R", "double"),
       no_sites = get_child(rep, "repl.no_sites", "double"),
       no_acts = get_child(rep, "repl.no_acts", "double"),
       type = get_child(rep, "repl.type", "double"),
       rev_type = get_child(rep, "repl.prev_type", "double")
  )
  acts <- get_child(rep, "activities")
  for(a in 1:xml_length(acts)){
    out[[paste0("act", a)]] <- xml_double(xml_child(acts, a))
  }
  return(out)
}))

rep_table
