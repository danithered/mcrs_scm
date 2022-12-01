library(shiny)
library(xml2)
library(DT)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(RRNA)
library(plotrix)

#Functions
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
  
  if(as.text) return( paste0("E[", paste(sort(acts), collapse = ","), "]") )
  return( as.expression(bquote(E[.(paste(sort(acts), collapse = ","))])) )
})
# Get data

rules <- readRDS("rules.RDS")

# Server
shinyServer(function(input, output) {
    setwd("/home/danielred/data/programs/mcrs_to_scm/OUT/A7retest.6_5/SAVE/")
    files <- grep(".xml", list.files(), value = T)
    
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
      selectizeInput("file",
                     label = "Which report would you like to see?",
                     #selectize = F,
                     options=list(maxOptions= length(files) ),
                     choices= files
      )
    })
    
    ## Reading in data
    observeEvent(input$file, {try({
      data <- read_xml(input$file)
      
      d <- xml_child(data, 1) # mooving to mcrscm
      
      time(get_child(d, "time", "int"))
      size(get_child(d, "sim.size", "int"))
      try(no_last_splits(get_child(d, "sim.no_last_splits", "int")))
      
      cells <- get_child(d, "cells")
      
      # get table
      tablev( do.call(rbind, lapply(1:xml_length(cells), function(no_cell){
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
      )
    }) })
    
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
        
        # refresh type colors
        #no_types = 2^no_acts_curr-1
        
        # refresh activity colors
        actcols( brewer.pal(no_acts_curr, "Set1") )
        
        # refresh replicator plot colorization
        pcols = list()
        for(col in actcols()) {
          pcols[[length(pcols)+1]] <- c("red", "coral")
        }
        col.pattern(pcols)
      }
      
      
    })
    
    ## Tables
    output$table <- renderDataTable( {
      tablev()[, -5]
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
      out
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
      rep_table()[whichone, c("seq") ]
    })
    
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
    output$hist_mfe <- renderPlot({
      
      #hist(rep_table()$mfe)
      ggplot(rep_table(), aes(x=mfe))+
        geom_histogram()
      
    })
    
    output$hist_Pfold <- renderPlot({
      ggplot(rep_table(), aes(x=Pfold))+
        geom_histogram()
    })
    
    output$hist_Pdeg <- renderPlot({
      ggplot(rep_table(), aes(x=Pdeg))+
        geom_histogram()
    })
    
    output$hist_R <- renderPlot({
      ggplot(rep_table(), aes(x=R))+
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
    
})
