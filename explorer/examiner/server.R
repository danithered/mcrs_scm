library(shiny)
library(xml2)
library(DT)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(RRNA)
library(plotrix)
library(dplyr)
library(Polychrome)
library(ssh)

#Functions
sshready_adr <- function(orig){
  unlist(lapply(strsplit(orig, ":"), function(x) paste(x, collapse=" -p ")))
}

dirto <- function(to, from=c(0,0)){
  
  return(atan2(to[2] - from[2], to[1] - from[1]) / pi)
}
rad2deg <- function(rad) {rad * 180 / pi}
pirad2deg <- function(rad) {rad * 180}
deg2rad <- function(deg) {deg * pi / 180}
deg2pirad <- function(deg) {deg / 180}
addToPlot <- function(x=0, y=0, coords, 
                      xspan=1, rot = NA, gap=0.05, 
                      main_con=NA, side_con=NA, 
                      add_letter=F, cex_letter=par("cex"), col_letter="black", col=NULL,
                      ...){
  par(bg="transparent")
  #rotate
  if(!is.na(rot)){
    midx <- mean(coords$x)
    midy <- mean(coords$y)
    dist <- sqrt(coords$x^2 + coords$y^2)
    angle <- atan2(coords$y, coords$x) + rot
    coords$y <- midy + sin(angle) * dist
    coords$x <- midx + cos(angle) * dist
    
  }
  
  #align
  coords$x <- (coords$x - min(coords$x))
  coords$y <- (coords$y - min(coords$y))
  
  #rescale
  scaling=1
  if(!is.na(xspan)){
    #compute scaling factor
    scaling = ifelse(max(coords$x) < max(coords$y),
                     xspan/(max(coords$y) - min(coords$y)),
                     xspan/(max(coords$x) - min(coords$x))
    )
    
    coords$x <- coords$x * scaling
    coords$y <- coords$y * scaling
  }
  
  #realign
  coords$x <- x + coords$x
  coords$y <- y + coords$y
  
  #connector lines
  ## MAIN_CON
  if(!is.list(main_con)) if(!is.na(main_con)){
    if(is.logical(main_con)){
      if(main_con) {
        main_con <- list(lwd=1, lty=1, col="black")
      }
    }
  }
  
  if(is.list(main_con)) with(main_con, {
    #points(coords$x, coords$y, type="c", col=col, lwd=lwd, lty=lty)
    segments(coords[1:(nrow(coords)-1),"x"], 
             coords[1:(nrow(coords)-1),"y"], 
             coords[2:nrow(coords),"x"], 
             coords[2:nrow(coords),"y"], 
             col=col, lwd=lwd, lty=lty)
  })
  
  ##SIDE_CON
  if(!is.list(side_con)) if(!is.na(side_con)){
    if(is.logical(side_con)) if(side_con) side_con <- list(lwd=1, lty=1, col="black")
  }  
  
  if(is.list(side_con)) with(side_con, {
    segs <- data.frame()
    
    for(base in 1:nrow(coords)){
      if(coords[base, "bound"] > 0 & base < coords[base, "bound"]){
        segs <- rbind(segs, data.frame(
          x0=coords[base, "x"], 
          y0=coords[base, "y"],
          x1=coords[coords[base, "bound"], "x"],
          y1=coords[coords[base, "bound"], "y"]  ))
      }
    }
    segments(segs$x0, segs$y0, segs$x1, segs$y1, lwd=lwd, col=col, lty=lty)
    
  })
  
  
  #points(coords$x, coords$y)
  # symbols(coords$x, coords$y, 
  #         circles=rep(scaling* (1/2-gap), nrow(coords)), 
  #         inches=F, add=T, ...)
  if(is.null(col)){
    cols <- rep("transparent", nrow(coords))
  } else {
    cols <- col
  }
  for(ci in 1:nrow(coords) ) draw.circle(coords$x[ci], coords$y[ci],
                                         radius = scaling* (1/2-gap), nrow(coords),
                                         col=cols[ci],
                                         ...)
  #add letter
  lett <- coords[add_letter, ]
  if(nrow(lett) > 0) {
    text(lett$x, lett$y, labels=lett$seq, cex=cex_letter, col=col_letter)
  }
  
  #return(coords)
}
plot_RNA <- function(coords, ...){
  plot.new()
  plot.window(asp=1, xlim=c(0,1), ylim=c(0,1), xpd=NA)
  addToPlot(0,0,coords, xspan=1, ...)
}
plot_acts <- function(a1, a2, col="grey"){
  length =max(length(a1), length(a2)) 
  if(length(a1) < length) a1 <- c(unlist(a1), rep(0, length-length(a1)))
  if(length(a2) < length) a2 <- c(unlist(a2), rep(0, length-length(a2)))
  plot.new()
  plot.window(asp=1, xlim=c(0,length), ylim=c(0,1), xpd=NA)
  rect(0:(length-1), 0, 1:length,1, 
       col = ifelse( rep(is.logical(a2), length), ifelse(a2, col, "white"), ifelse(a2 == 0, "white", alpha(col, a1/max(a1)))))
  rect(0:(length-1), 1.5, 1:length,2.5, 
       col = ifelse( rep(is.logical(a1), length), ifelse(a1, col, "white"), ifelse(a1 == 0, "white", alpha(col, a1/max(a1)))))
  if(!is.logical(a1)) text(1:length-0.5, 2, labels = round(a1, 2) )
  if(!is.logical(a2)) text(1:length-0.5, 0.5, labels = round(a2, 2) )
  text(mean(c(0, length) ), 2.75, "Current strand")
  text(mean(c(0, length) ), 1.25, "Complementary strand")
}

find_activity <- function(repl, coords, rules){
  startofPattern <- list()
  for(rule in rules){
    where <- gregexpr(rule$str, repl$str, fixed=T)
    if(where[[1]] > -1){ #it has this activity
      for(w in where){ #check every possible location
        #check subrules
        found = c()
        for(sr in rule$subrules){
          at = w + as.numeric(sr$locus)-1
          if(substr(repl$seq, at, at) == sr$base) { 
            found = c(found, at)
          }
        } # going tru subrules
        if(length(found)>0){ #everything is ok
          startofPattern[[length(startofPattern)+1]] <- list(activity = rule$activity, 
                                                             start=as.numeric(w), 
                                                             length= nchar(rule$str), 
                                                             matches=found
          )
        }
      } # checking possible locations
    }
  }
  if(length(startofPattern) == 0) return(NA)
  return(startofPattern)
}
basecol <- function(str, n, col){
  sapply(n, function(n=n, str, col){
    if(length(col) == 1){
      return(col)
    }
    if(length(col) == 2){
      base <- substr(str, n, n)
      if(base == "."){
        return(col[1])
      }
      if(base %in% c("(", ")")){
        return(col[2])
      }
    }
    return(NA)
  }, str=str, col=col)
}
masking <- function(str, n, col){
  mask <- rep(NA, nchar(str))
  mask[n] <- basecol(str, n, col)
  return(mask)
}
mask_overlap <- function(mask1, mask2){
  vals2 = !is.na(mask2)
  overlap = !is.na(mask1) & vals2
  mask1[vals2] <- mask2[vals2]
  mask1[overlap] <- "orange"
  return(mask1)
}
make.colormask <- function(str, patterns=list(), col=NA, col.pattern=list(), col.base=c()){
  #browser()
  length = nchar(str)
  basemask = basecol(str, 1:length, col)
  if(!is.na(patterns)[1]) for(pattern in patterns){
    if(length(col.pattern) < pattern$activity) mask <- rep(NA, length)
    else mask <- masking(str, seq(pattern$start, length.out=pattern$length), col.pattern[[pattern$activity]])
    if(length(col.base) >= pattern$activity) mask[pattern$matches] <- col.base[pattern$activity]
    basemask <- mask_overlap(basemask, mask)
  }
  return(basemask)
}


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

type2noA <- Vectorize(function(x){
  if(x==0) return("A 0")
  if(x==-1) return("empty")
  
  m <- ceiling(log(x,2))
  acts <- c()
  
  for(i in m:0){
    if(x %/% 2^i != 0){
      acts <- c(acts, i)
      x <- x-2^i
      if(x == 0) break
    }
  }
  return( paste("A",length(acts)) )
  
})

type2A <- Vectorize(function(x){
  if(x==0) return(0)
  if(x==-1) return(NA)
  
  m <- ceiling(log(x,2))
  acts <- c()
  
  for(i in m:0){
    if(x %/% 2^i != 0){
      acts <- c(acts, i)
      x <- x-2^i
      if(x == 0) break
    }
  }
  return( length(acts) )
  
})
type2vec <- (function(x, no=1){
  if(x==0) return(rep(F, no))
  if(x==-1) return(NA)
  
  m <- ceiling(log(x,2))
  acts <- c()
  
  for(i in m:0){
    if(x %/% 2^i != 0){
      acts <- c(acts, i)
      x <- x-2^i
      if(x == 0) break
    }
  }
  
  acts <- acts+1
  out <- rep(F, max(no, acts))
  out[acts] <- T
  
  return(out)
})

enzN <- Vectorize(function(x, as.text=F){
  if(x==0) return("parazite")
  if(x==-1) return("empty")
  
  m <- ceiling(log(x,2))
  acts <- c()
  
  for(i in m:0){
    if(x %/% 2^i != 0){
      acts <- c(acts, i)
      x <- x-2^i
      if(x == 0) break
    }
  }
  
  if(as.text) return( paste0("E[", paste(sort(acts+1), collapse = ","), "]") )
  return( as.expression(bquote(E[.(paste(sort(acts+1), collapse = ","))])) )
})

mergepath <- function(...){
  x <- as.character(unlist(list(...)))
  x <- x[nchar(x) > 0]
  if(length(x) == 0) return(NA)
  
  wrong <- c(substr(x[1:(length(x)-1) ], nchar(x), nchar(x)) != "/", F)
  x[wrong] <- paste0(x[wrong], "/")
  
  wrong <- c(F, substr(x[2:length(x)], 1, 1) == "/")
  x[wrong] <- substr(x[wrong], 2, nchar(x[wrong]))
  
  return( paste(x, collapse = "") )
}

get_filelist <- Vectorize(function(path="~/", ssh=NA, ssh_key="~/.ssh/id_rsa"){
  if(nchar(path) > 0 & substr(path, nchar(path), nchar(path)) != "/" ) path=paste0(path, "/")
  
  if(is.na(ssh)){ #local dir
    if( dir.exists( path ) ){
      out <- system(paste("find", path, "-maxdepth 1 -type f"), intern=T)
      if( length(out) == 0 ) return( NA) 
      else return( out )
    } else {
      return(NA)
    }
  } else try({ 
    #it is in shh
    #connecting to ssh
    ssh_con <- ssh_connect(ssh, keyfile = ssh_key)
    if(!ssh_session_info(ssh_con)$connected){
      warning("Could not establish shh connection")
      return(NA)
    }
    
    #check if dir exists
    if(!as.logical(capture.output(ssh_exec_wait(ssh_con, paste0("if [ -d ", path, " ]; then echo TRUE; else echo FALSE; fi")  ))[1])){
      return(NA)
    }
    
    #get subdirs
    out <- capture.output(ssh_exec_wait(ssh_con, paste("find", path, "-maxdepth 1 -type f") ))
    if(length(out) > 1) out <- out[1:(length(out)-1)]
    else out <- NA
    
    #disconnect
    ssh_disconnect(ssh_con)
    return(out)
    
  })
})

get_file <- function(target, path="~/", ssh=NA, ssh_key="~/.ssh/id_rsa", fast=T, to=NA){
  if(is.na(ssh)){
    file = paste(path, 
                 target, 
                 sep=ifelse(nchar(path) > 0 & substr(path, nchar(path), nchar(path)) != "/", "/", ""))
    if( all(file.exists(file)) ) return(file)
    else return(NA)
  } else try({ 
    #it is in shh
    
    #check output dir
    if(is.na(to)) {
      to = tempdir()
    } else {
      if(!dir.exists(to)){
        if(!dir.create(to)) {
          warning(paste("Cannot create dir", to))
          return(NA)
        }
      }
    }
    
    #download
    tofile = paste(to, filename(target), sep=ifelse(nchar(to) > 0 & substr(to, nchar(to), nchar(to)) != "/", "/", ""))
    if(fast){
      system(paste("ssh", sshready_adr(ssh), '"tar -C', path, '-zc', target, '" | tar zx -C', to))
      if(all(file.exists( tofile ))){
        return( tofile )
      } else {
        warning(paste("Could not download file tru shh", target))
        return(NA)
      }
    } else {
      #connecting to ssh
      ssh_con <- ssh_connect(ssh, keyfile = ssh_key)
      if(!ssh_session_info(ssh_con)$connected){
        warning("Could not establish shh connection")
        return(NA)
      }
      #download
      scp_download(ssh_con, 
                   files= paste(path, target, sep=ifelse(nchar(path) > 0 & substr(path, nchar(path), nchar(path)) != "/", "/", "")),
                   to=to
      )
      #disconnect
      ssh_disconnect(ssh_con)
      return( tofile )
    }
  })
}

filename <- function(x){
  sapply(strsplit(x, "/"), function(x) x[length(x)])
}

path <- function(x){
  sapply(strsplit(x, "/"), function(x) ifelse(length(x) > 1, paste(x[1:(length(x)-1)], collapse = "/"), x) )
}

# app specific functions
get_my_data <- function(file, path="", ssh=NA, ssh_key= "~/.ssh/id_rsa"){
  out <- list()
  
  try({
    f <- get_file(file, path, ssh=ssh, ssh_key=ssh_key, fast=T, to=NA)
  
    data <- read_xml(f)
    
    d <- xml_child(data, 1) # mooving to mcrscm
    
    out$time = get_child(d, "time", "int")
    out$size = get_child(d, "sim.size", "int")
    try({
      out$no_last_splits = get_child(d, "sim.no_last_splits", "int")
      })
    
    cells <- get_child(d, "cells")
    
    # get table
    out$table = do.call(rbind, lapply(1:xml_length(cells), function(no_cell){
      cell <- xml_child(cells, no_cell)
      
      out <- as.data.frame(list( alive = get_child(cell, "cell.alive", "logical"),
                                 leftover = get_child(cell, "cell.leftover", "double"),
                                 M = get_child(cell, "metabolism", "double"),
                                 no_reps = get_child(cell, "cell.reps") |> get_child("count", "int")#,
                                 #reps = get_child(cell, "cell.reps")
      ))
      out$reps <- list(get_child(cell, "cell.reps"))
      out
    }))
    
    if(!is.na(ssh)) file.remove(f)
  }) # try
  
  return(out)
}

get_my_xmls <- function(path, ssh=NA, ssh_key="~/.ssh/id_rsa"){
  grep(".xml", get_filelist(path= mergepath(path, "SAVE/"), ssh=ifelse(nchar(ssh) == 0, NA, ssh), ssh_key=ssh_key), value=T)
}

# Get data

rules <- readRDS("rules.RDS")

# Server
shinyServer(function(input, output, session) {
    setwd("/home/danielred/data/programs/mcrs_to_scm/OUT/A7retest.6_5/SAVE/")
  
    # inic params
    params <- reactiveValues()
    params$cache.path <- "report_cache/"
    params$dir <- "/home/danielred/data/programs/mcrs_to_scm/OUT/A7retest.6_5/"
    params$ssh <- NA
    params$ssh_key <- "~/.ssh/id_rsa"
    params$force <- FALSE
    
    # reads params from URL
    observe({
      query <- parseQueryString(session$clientData$url_search)
      if (length(query) > 0) {
        ws <- names(query)[names(query) %in% names(params)]
        for(w in ws) params[[w]] <- ifelse(nchar(query[[w]]) > 0, query[[w]], NA)
        files(get_my_xmls( params$dir, params$ssh, params$ssh_key))
      }
    })
  
    # inic other vals
    files <- reactiveVal(get_my_xmls( isolate(params$dir), isolate(params$ssh), isolate(params$ssh_key))) 
    tablev <- reactiveVal()
    rep_table <- reactiveVal()
    size <- reactiveVal()
    time <- reactiveVal()
    no_last_splits <- reactiveVal()
    typecols <- reactiveVal()
    actcols <- reactiveVal()
    no_acts <- reactiveVal(0)
    col.pattern <- reactiveVal()
      
    ## UI
    output$reports <- renderUI({
      files()
      fs <- isolate(files())
      names(fs) <- filename(fs)
      nums <- substr(names(fs), 1, nchar(names(fs))-4)
      if(all(!is.na(as.numeric(nums)))){ # if all is numeric
        #names(files) <- sapply(strsplit(names(files), ".", fixed=T), function(x) x[1])
        names(fs) <- as.character(nums)
        selectizeInput("file",
                       label = "Which report would you like to see?",
                       #selectize = F,
                       options=list(maxOptions= length(fs)),
                       choices= fs[order(as.numeric(nums))]
        )
      } else { # there are non numeric ones
        selectizeInput("file",
                       label = "Which report would you like to see?",
                       #selectize = F,
                       options=list(maxOptions= length(fs )),
                       choices= fs
        )
      }
    })
    
    popup <- function(pdir, pssh, pssh_key) modalDialog(
      textInput("path", "Path:", value=pdir),
      textInput("ssh", "SSH:", value=pssh),
      textInput("ssh_key", "PAth to SSH key:", value = pssh_key),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("submit_path", "OK")
      )
    )
    
    # Reading in data
    observeEvent(input$pop, {
      showModal(popup(params$dir, params$ssh, params$ssh_key))
    })
    
    observeEvent(input$submit_path, {
      params$dir <- input$path
      params$ssh <- ifelse(nchar(input$ssh) > 0, input$ssh, NA)
      params$ssh_key <- input$ssh_key
      
      files(get_my_xmls( params$dir, params$ssh, params$ssh_key))
      
      removeModal()
    })
    
    observeEvent(input$file, {
      data <- get_my_data(filename(input$file), path=path(input$file), ssh=params$ssh, ssh_key = params$ssh_key)
      
      if(length(data) > 0){
        no_last_splits(data$no_last_splits)
        time(data$time)
        size(data$size)
        tablev(data$table)
      }
    }) # observeEvent
    
    observeEvent(input$table_rows_selected, {
      # examine reps
      reps = tablev()[as.numeric(input$table_rows_selected),]$reps[[1]]
      
      rep_table (do.call(rbind, lapply(which(children(reps) == "item"), function(rnum){
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
        return(as.data.frame(out))
      })))
      
      # create colors
      no_acts_curr <- sum(startsWith(colnames(rep_table()), "act"))
      if(no_acts_curr != no_acts()){ # number of activities changed
        # refresh no_acts
        no_acts(no_acts_curr)
        
        # refresh activity colors
        actcols( brewer.pal(no_acts_curr, "Set1") )
        
        # refresh replicator plot colorization
        pcols = list()
        for(col in actcols()) {
          pcols[[length(pcols)+1]] <- c("red", "coral")
        }
        col.pattern(pcols)
        
        # refresh type colors
        no_types = 2^no_acts_curr-1
        no_types
        set.seed(5464)
        newcols <- c(NA, col2rgb("black"), createPalette(no_types, actcols()))
        names(newcols) <- enzN(-1:no_types, as.text = T)
        newcols[!is.na(type2A(-1:no_types)) & type2A(-1:no_types)==1] <- actcols()
        typecols(newcols)
      }
      
      
    })
    
    ## Tables
    output$table <- renderDataTable( {
      tablev()[, -5] |> mutate_if(is.numeric, round, digits=2)
    }, 
    selection =list(mode = 'single', selected = 1, target = 'row', selectable = T), 
    options = list(
      order = list(list(1, "desc"), list(3, "desc")),
      scrollY = 200,
      scroller = TRUE,
      deferRender = TRUE
    ),
    extensions = "Scroller"
    )
    
    output$reps <- renderDataTable({
      out <- cbind(
        type = enzN(rep_table()$type, as.text=T),
        reverse = enzN(rep_table()$rev_type, as.text = T),
        length = nchar(rep_table()$seq),
        rep_table()[,c("mfe", "Pfold", "Pdeg", "R", "no_sites", "no_acts")],
        rep_table()[, startsWith(colnames(rep_table()), "act")]
      )
      out |> mutate_if(is.numeric, round, digits=2)
    }, 
    selection =list(mode = 'single', selected = 1, target = 'row', selectable = T), 
    options = list(
      scrollY = 200,
      scroller = TRUE,
      deferRender = TRUE
    ),
    extensions = "Scroller"
    )
    
    output$rep_props <- renderTable({
      whichone= input$reps_rows_selected
      t(as.data.frame(c(length=nchar(rep_table()[whichone, "seq"]), rep_table()[whichone, c("mfe", "Pfold", "Pdeg", "R") ]))) |> round(digits=2)
    }, rownames = T, colnames = F)
    
    ## Text outputs
    output$state_data <- renderText(paste0("Sample taken in generation: ", 
                                           time(), 
                                           "\nNumber of cells: ", size(), 
                                           "\nNumber of splitting events in last generation: ",
                                           no_last_splits() )) 
    output$seq <- renderText({
      whichone= input$reps_rows_selected
      rep_table()[whichone, c("seq") ]
    })
    output$str <- renderText({
      whichone= input$reps_rows_selected
      rep_table()[whichone, c("str") ]
    })
    
    ## Plots
    # For state
    output$plot <- renderPlotly( {
      s = input$table_rows_selected
      #plot(tablev()$M, tablev()$no_reps)
      #points(tablev()[s,"M"], tablev()[s, "no_reps"], col="red")
      
      sizes = rep(1,nrow(tablev()))
      sizes[s] <- 4
      ggplotly(ggplot(tablev(), aes(x=M, y= no_reps, color=alive))+
        geom_point(size=sizes )+
        labs(x="Metabolism", y="Number of replicators")+
        theme(legend.pos="none")
      )
    })
    
    # For cell
    output$legend <- renderPlot({
      plot.new()
    })
    
    output$hist_mfe <- renderPlot({
      
      #hist(rep_table()$mfe)
      ggplot(rep_table(), aes(x=mfe, fill=as.factor(type)))+
        geom_histogram()
      
    })
    
    output$hist_Pfold <- renderPlot({
      ggplot(rep_table(), aes(x=Pfold))+
        geom_histogram()+
        scale_fill_manual(values=typecols())+
        theme(legend.position = "none")
    })
    
    output$hist_Pdeg <- renderPlot({
      ggplot(rep_table(), aes(x=Pdeg))+
        geom_histogram()
    })
    
    output$hist_R <- renderPlot({
      ggplot(rep_table(), aes(R, fill=no_acts))+
        geom_histogram()
    })
    
    output$hist_no_sites <- renderPlot({
      ggplot(rep_table(), aes(x=no_sites))+
        geom_histogram()
    })
    
    output$hist_no_acts <- renderPlot({
      ggplot(rep_table(), aes(x=no_acts))+
        geom_histogram()
    })
    
    output$nice <- renderPlot({
      types <- unique(rep_table()$type)
      p <- list(par_noEA=7)
      
      # calculate table
      van <- rep_table()$seq != "N"
      odf <- data.frame(orig=rep_table()$type[van], rev=rep_table()$rev_type[van])
      odf$orig <- factor(odf$orig, levels = 0:(2^as.numeric(p$par_noEA)-1))
      odf$rev <- factor(odf$rev, levels = 0:(2^as.numeric(p$par_noEA)-1))
      pairs <- table(odf)
      
      tv <- as.numeric(colnames(pairs))
      
      kell <- 0:(2^as.numeric(p$par_noEA)-1)

      #reorder by numbers of activities
      kell2 <- kell[order(type2noA(kell))]
      pairs <- pairs[as.character(kell2),]
      pairs <- pairs[,as.character(kell2)]
      
      keep <- apply(pairs, 1, function(x) sum(x) > 0 ) | apply(pairs, 2, function(x) sum(x) > 0 )
      kell <- kell[keep]
      pairs <- pairs[keep, keep]
    
      image(1:sum(keep), 1:sum(keep),  
            pairs, 
            xaxt="n", yaxt="n" ,
            xlab="", ylab="",
            col= heat.colors(100, rev = TRUE)[5:100]
      )
      axis(1, at= 1:sum(keep), labels=enzN(kell), las=2)
      axis(2, at= 1:sum(keep), labels=enzN(kell), las=1)
      text(rep(1:sum(keep), length(1:sum(keep))), rep(1:sum(keep), each=length(1:sum(keep))), c(pairs), cex=0.5 )
      
    })

    # For Replicatro
    output$replicator <- renderPlot({
      par(mar=c(0,0,0,0))
      whichone= input$reps_rows_selected
      
      repl <- rep_table()[whichone, c("seq", "str") ]
      
      sink("NUL") # they made ct2coord to print everything
      coords= ct2coord( makeCt( repl$str, repl$seq) )
      sink()
      
      startofPatterns <- find_activity(repl, coords, rules)
      
      #compute 2D structure and colorisation
      colormask <- make.colormask(repl$str, 
                                  patterns =startofPatterns, 
                                  col=NA, 
                                  col.pattern = col.pattern(), 
                                  col.base = actcols()
      )
      
      plot_RNA(coords,
               border="lightblue",
               #bases
               add_letter = T,
               cex_letter = 0.6,
               col_letter = "black",
               #fill
               col=colormask,
               #rotate it
               #rot=pi/2,
               #connecting lines
               main_con = list(lwd=1, col="darkgrey", lty=1),
               side_con = list(lwd=0.5, col="purple", lty=2)
      )
      
    }) # render replicator
    
    output$acts <- renderPlot({
      par(mar=c(0,0,0,0))
      whichone= input$reps_rows_selected
      
      acts = rep_table()[whichone, ] |> select( starts_with("act"))
      compl_acts = type2vec(rep_table()[whichone, "rev_type"], no_acts())
      plot_acts(acts, compl_acts, col=actcols())
    })
    
})
