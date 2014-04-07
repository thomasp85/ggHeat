ggHeat <- function(data, colfactor, rowfactor, cGroup, rGroup, cOrder, rOrder, hideticks='none', dataname='', method='ward', distance='euclidean', title='', outputPlot=FALSE, legend=FALSE, nPages=1, suppress.leg, ...){
    if(length(method) == 1){
        method <- rep(method, 2)
    } else {}
    if(length(distance) == 1){
        distance <- rep(distance, 2)
    } else {}
    thetimer <- list()
    
    # Save dimnames for later and put in unique numbers
    rawDimNames <- list(col=colnames(data), row=rownames(data))
    colnames(data) <- 1:ncol(data)
    rownames(data) <- 1:nrow(data)
    
    # Create clustering
    
    cat('Creating clustering of data...\n')
    flush.console()
    
    clustwithinFac <- function(x, sort){
        sort <- sort[x$sort]
        ord <- order(order(unique(sort)))[match(sort, unique(sort))]
        x$order <- as.numeric(paste(x$order, sprintf(paste('%0', nchar(length(ord)), 'd', sep=''), ord), sep=ifelse(any(grepl('.', x$order, fixed=TRUE)), '', '.')))
        x
    }
    clustwithinAuto <- function(x, data){
        if(nrow(x) > 1){
            data <- data[, x$sort]
            ord <- order(hclust(dist(t(data), method=distance[2]), method=method[2])$order)
            x$order <- as.numeric(paste(x$order, sprintf(paste('%0', nchar(length(ord)), 'd', sep=''), ord), sep=ifelse(any(grepl('.', x$order, fixed=TRUE)), '', '.')))
        } else {}
        x
    }
    
    ## Rows
    start <- proc.time()
    if(missing(rOrder)){
        rows <- hclust(dist(data, method=distance[1]), method=method[1])
        if(!missing(rGroup)){
            if(rGroup < 1){
                height <- max(rows$height)*(1-rGroup)
                rGroup <- cutree(rows, h=height)
            } else {
                rGroup <- cutree(rows, k=rGroup)
            }
        } else {}
        rows <- dendro_data(rows)
        if(!missing(rGroup)){
            rows <- denCut(rows, rGroup)
        } else {}
        rowlab <- as.character(rows$labels$label)
        rOrdertype <- 'none'
    } else {
        if(!all(rOrder %in% names(rowfactor))){
            stop('Factor to sort by not part of supplied data frame')
        } else {}
        rows <- rowfactor[,rOrder[1]]
        if(is.numeric(rows)){
            rOrdertype <- 'seq'
            rOrdering <- order(order(unique(rows)))[match(rows, unique(rows))]
            rOrdering <- data.frame(factor=rows, sort=1:length(rOrdering), order=rOrdering)
            rOrdertmp <- rOrder[-1]
            while(length(rOrdertmp) > 0){
                rOrderingsub <- rowfactor[, rOrdertmp[1]]
                rOrdering <- ddply(rOrdering, .(order), function(x, sort) clustwithinFac(x, sort), sort=rOrderingsub)
                rOrdering <- rOrdering[order(rOrdering$sort), ]
                rOrdertmp <- rOrdertmp[-1]
            }
            rOrdering <- ddply(rOrdering, .(order), function(x, data) clustwithinAuto(x, t(data)), data=data)
            rOrdering <- rOrdering[order(rOrdering$sort), ]
            rOrdering <- order(rOrdering$order)
            rowlab <- as.character(rownames(data)[rOrdering])
        } else {
            rOrdertype <- 'qual'
            rOrdering <- rows
            rowfactor <- rowfactor[, which(names(rowfactor) != rOrder[1]), drop=FALSE]
            rOrderdat <- ddply(data.frame(factor=rows, data), .(factor), function(x) apply(x[,-1], 2, function(x) mean(x, na.rm=T)))
            rOrderlev <- rOrderdat[,1]
            rOrderdat <- rOrderdat[,-1]
            rows <- hclust(dist(rOrderdat, method=distance[1]), method=method[1])
            rOrderlev <- data.frame(factor=rOrderlev[rows$order], order=1:length(rOrderlev))
            rOrdering <- merge(data.frame(factor=rOrdering, sort=1:length(rOrdering)), rOrderlev, sort=FALSE)
            rOrdering <- rOrdering[order(rOrdering$sort), ]
            rOrdering$tmp <- rOrdering$order
            rOrdertmp <- rOrder[-1]
            while(length(rOrdertmp) > 0){
                rOrderingsub <- rowfactor[, rOrdertmp[1]]
                rOrdering <- ddply(rOrdering, .(order), function(x, sort) clustwithinFac(x, sort), sort=rOrderingsub)
                rOrdering <- rOrdering[order(rOrdering$sort), ]
                rOrdertmp <- rOrdertmp[-1]
            }
            rOrdering <- ddply(rOrdering, .(order), function(x, data) clustwithinAuto(x, t(data)), data=data)
            rOrdering <- rOrdering[order(rOrdering$sort), ]
            rOrderlev <- merge(rOrderlev, ddply(rOrderlev, .(order), function(x) range(which(rOrdering[order(rOrdering$order), 'tmp'] == x[1, 'order']))))
            names(rOrderlev)[3:4] <- c('min', 'max')
            rOrdering <- order(rOrdering$order)
            rowlab <- as.character(rownames(data)[rOrdering])
            rows <- dendro_data(rows)
            rows$labels$label <- rOrderlev$factor
        }
    }
    thetimer$Cluster.rows <- (proc.time() - start)[3]
    start <- proc.time()
    
    ## Columns
    if(missing(cOrder)){
        cols <- hclust(dist(t(data), method=distance[2]), method=method[2])
        if(!missing(cGroup)){
            if(cGroup < 1){
                height <- max(cols$height)*(1-cGroup)
                cGroup <- cutree(cols, h=height)
            } else {
                cGroup <- cutree(cols, k=cGroup)
            }
        } else {}
        cols <- dendro_data(cols)
        if(!missing(cGroup)){
            cols <- denCut(cols, cGroup)
        } else {}
        collab <- as.character(cols$labels$label)
        cOrdertype <- 'none'
    } else {
        if(!all(cOrder %in% names(colfactor))){
            stop('Factor to sort by not part of supplied data frame')
        } else {}
        cols <- colfactor[,cOrder[1]]
        if(is.numeric(cols)){
            cOrdertype <- 'seq'
            cOrdering <- order(order(unique(cols)))[match(cols, unique(cols))]
            cOrdering <- data.frame(factor=cols, sort=1:length(cOrdering), order=cOrdering)
            cOrdertmp <- cOrder[-1]
            while(length(cOrdertmp) > 0){
                cOrderingsub <- colfactor[, cOrdertmp[1]]
                cOrdering <- ddply(cOrdering, .(order), function(x, sort) clustwithinFac(x, sort), sort=cOrderingsub)
                cOrdering <- cOrdering[order(cOrdering$sort), ]
                cOrdertmp <- cOrdertmp[-1]
            }
            cOrdering <- ddply(cOrdering, .(order), function(x, data) clustwithinAuto(x, data), data=data)
            cOrdering <- cOrdering[order(cOrdering$sort), ]
            cOrdering <- order(cOrdering$order)
            collab <- as.character(colnames(data)[cOrdering])
        } else {
            cOrdertype <- 'qual'
            cOrdering <- cols
            colfactor <- colfactor[, which(names(colfactor) != cOrder[1]), drop=FALSE]
            cOrderdat <- ddply(data.frame(factor=cols, t(data)), .(factor), function(x) apply(x[,-1], 2, function(x) mean(x, na.rm=T)))
            cOrderlev <- cOrderdat[,1]
            cOrderdat <- cOrderdat[,-1]
            cols <- hclust(dist(cOrderdat, method=distance[2]), method=method[2])
            cOrderlev <- data.frame(factor=cOrderlev[cols$order], order=1:length(cOrderlev))
            cOrdering <- merge(data.frame(factor=cOrdering, sort=1:length(cOrdering)), cOrderlev, sort=FALSE)
            cOrdering <- cOrdering[order(cOrdering$sort), ]
            cOrdering$tmp <- cOrdering$order
            cOrdertmp <- cOrder[-1]
            while(length(cOrdertmp) > 0){
                cOrderingsub <- colfactor[, cOrdertmp[1]]
                cOrdering <- ddply(cOrdering, .(order), function(x, sort) clustwithinFac(x, sort), sort=cOrderingsub)
                cOrdering <- cOrdering[order(cOrdering$sort), ]
                cOrdertmp <- cOrdertmp[-1]
            }
            cOrdering <- ddply(cOrdering, .(order), function(x, data) clustwithinAuto(x, data), data=data)
            cOrdering <- cOrdering[order(cOrdering$sort), ]
            cOrderlev <- merge(cOrderlev, ddply(cOrderlev, .(order), function(x) range(which(cOrdering[order(cOrdering$order), 'tmp'] == x[1, 'order']))))
            names(cOrderlev)[3:4] <- c('min', 'max')
            cOrdering <- order(cOrdering$order)
            collab <- as.character(colnames(data)[cOrdering])
            cols <- dendro_data(cols)
            cols$labels$label <- cOrderlev$factor
        }
    }
    thetimer$Cluster.cols <- (proc.time() - start)[3]
    start <- proc.time()
    
    # Transform factor info
    
    cat('Transforming factor information...\n')
    flush.console()
    
    ## Columns
    if(!missing(colfactor)){
        if(length(colfactor) != 0){
            colfacname <- names(colfactor)
            colfac <- list()
            for(i in 1:length(colfacname)){
                if(is.numeric(colfactor[,i])){
                    type <- 'seq'
                    colfactor1 <- colfactor[match(collab, colnames(data)), i]
                    colfactor1 <- binNum(colfactor1, 9)
                    names(colfactor1) <- c('value', 'colour')
                    colfactor1 <- data.frame(x=1:nrow(colfactor1), y=-1000*i, colfactor1, gridrow=1+i, gridcol=1, type=type)
                    colfactor1 <- split(colfactor1, colfactor1$colour)
                } else {
                    type <- 'qual'
                    colfactor1 <- as.character(colfactor[match(collab, colnames(data)), i])
                    colfactor1 <- data.frame(x=1:length(colfactor1), y=-1000*i, colour=colfactor1, gridrow=1+i, gridcol=1, type=type)
                    colfactor1 <- split(colfactor1, colfactor1$colour)
                }
                colfac[[i]] <- colfactor1
            }
        } else colfacname <- c()
    } else {
        colfacname <- c()
    }
    thetimer$Format.col.factor <- (proc.time() - start)[3]
    start <- proc.time()
    
    ## Rows
    if(!missing(rowfactor)){
        if(length(rowfactor) != 0){
            rowfacname <- names(rowfactor)
            rowfac <- list()
            for(i in 1:length(rowfacname)){
                if(is.numeric(rowfactor[,i])){
                    type <- 'seq'
                    rowfactor1 <- rowfactor[match(rowlab, rownames(data)), i]
                    rowfactor1 <- binNum(rowfactor1, 9)
                    names(rowfactor1) <- c('value', 'colour')
                    rowfactor1 <- data.frame(y=1:nrow(rowfactor1), x=-1000*i, rowfactor1, gridcol=1+i, gridrow=2+length(colfacname), type=type)
                    rowfactor1 <- split(rowfactor1, rowfactor1$colour)
                } else {
                    type <- 'qual'
                    rowfactor1 <- as.character(rowfactor[match(rowlab, rownames(data)), i])
                    rowfactor1 <- data.frame(y=1:length(rowfactor1), x=-1000*i, colour=rowfactor1, gridcol=1+i, gridrow=2+length(colfacname), type=type)
                    rowfactor1 <- split(rowfactor1, rowfactor1$colour)
                }
                rowfac[[i]] <- rowfactor1
            }
        } else rowfacname <- c()
    } else {
        rowfacname <- c()
    }
    thetimer$Format.row.factor <- (proc.time() - start)[3]
    start <- proc.time()
    
    ## Create color scales
    
    ### Factor colour
    seqname <- c('Greens', 'Oranges', 'Purples', 'Greys', 'Blues', 'Reds')
    seqno <- 1
    two <- 1
    three <- 1
    four <- 1
    
    ### for columns
    if(!missing(colfactor) & length(colfacname)!=0){
        for(i in 1:length(colfac)){
            if(colfac[[i]][[1]]$type[1] == 'seq'){
                bpal <- seqname[seqno]
                if(length(colfac[[i]]) == 2){
                    colour <- brewer.pal(3, bpal)[c(1,3)]
                } else {
                    colour <- brewer.pal(length(colfac[[i]]), bpal)
                }
                if(seqno == 6){
                    seqno <- 1
                } else {
                    seqno <- seqno+1
                }
            } else {
                if(length(colfac[[i]]) == 2){			  
                    if(two <= 4){
                        colour <- brewer.pal(8, 'Dark2')[c(1,2)+2*(two-1)]
                        bpal <- 'Dark2'
                    } else if (two <= 8){
                        colour <- brewer.pal(8, 'Accent')[c(1,2)+2*(two-5)]
                        bpal <- 'Accent'
                    } else {}
                    two <- two + 1
                    if(two == 9){
                        two <- 1
                    } else {}
                } else if(length(colfac[[i]]) == 3){			  
                    if(three <= 3){
                        colour <- brewer.pal(9, 'Set1')[1:3+3*(three-1)]
                        bpal <- 'Set1'
                    } else if (three <= 7){
                        colour <- brewer.pal(12, 'Set3')[1:3+3*(three-4)]
                        bpal <- 'Set3'
                    } else {}
                    three <- three + 1
                    if(three == 8){
                        three <- 1
                    } else {}
                } else if(length(colfac[[i]]) == 4){			  
                    if(four <= 2){
                        colour <- brewer.pal(8, 'Set2')[1:4+4*(four-1)]
                        bpal <- 'Set2'
                    } else if (four <= 5){
                        colour <- brewer.pal(12, 'Paired')[1:4+4*(four-3)]
                        bpal <- 'Paired'
                    } else {}
                    four <- four + 1
                    if(three == 6){
                        three <- 1
                    } else {}
                } else if(length(colfac[[i]]) > 12){
                    colour <- rainbow(length(colfac[[i]]))
                    bpal <- 'Rainbow'
                } else {
                    colour <- brewer.pal(length(colfac[[i]]), 'Paired')
                    bpal <- 'Paired'
                }
            }
            for(j in 1:length(colfac[[i]])){
                colfac[[i]][[j]]$colour <- colour[j]
                colfac[[i]][[j]]$Palet <- bpal
            }
        }
    } else {}
    thetimer$Colour.col.factor <- (proc.time() - start)[3]
    start <- proc.time()
    
    ### for Rows
    if(!missing(rowfactor) & length(rowfacname)!=0){
        for(i in 1:length(rowfac)){
            if(rowfac[[i]][[1]]$type[1] == 'seq'){
                bpal <- seqname[seqno]
                if(length(rowfac[[i]]) == 2){
                    colour <- brewer.pal(3, seqname[seqno])[c(1,3)]
                } else {
                    colour <- brewer.pal(length(rowfac[[i]]), seqname[seqno])
                }
                if(seqno == 6){
                    seqno <- 1
                } else {
                    seqno <- seqno+1
                }
            } else {
                if(length(rowfac[[i]]) == 2){			  
                    if(two <= 4){
                        colour <- brewer.pal(8, 'Dark2')[c(1,2)+2*(two-1)]
                        bpal <- 'Dark2'
                    } else if (two <= 8){
                        colour <- brewer.pal(8, 'Accent')[c(1,2)+2*(two-5)]
                        bpal <- 'Accent'
                    } else {}
                    two <- two + 1
                    if(two == 9){
                        two <- 1
                    } else {}
                } else if(length(rowfac[[i]]) == 3){			  
                    if(three <= 3){
                        colour <- brewer.pal(9, 'Set1')[1:3+3*(three-1)]
                        bpal <- 'Set1'
                    } else if (three <= 7){
                        colour <- brewer.pal(12, 'Set3')[1:3+3*(three-4)]
                        bpal <- 'Set3'
                    } else {}
                    three <- three + 1
                    if(three == 8){
                        three <- 1
                    } else {}
                } else if(length(rowfac[[i]]) == 4){			  
                    if(four <= 2){
                        colour <- brewer.pal(8, 'Set2')[1:4+4*(four-1)]
                        bpal <- 'Set2'
                    } else if (four <= 5){
                        colour <- brewer.pal(12, 'Paired')[1:4+4*(four-3)]
                        bpal <- 'Paired'
                    } else {}
                    four <- four + 1
                    if(three == 6){
                        three <- 1
                    } else {}
                } else if(length(rowfac[[i]]) > 12){
                    colour <- rainbow(length(rowfac[[i]]))
                    bpal <- 'Rainbow'
                } else {
                    colour <- brewer.pal(length(rowfac[[i]]), 'Paired')
                    bpal <- 'Paired'
                }
            }
            for(j in 1:length(rowfac[[i]])){
                rowfac[[i]][[j]]$colour <- colour[j]
                rowfac[[i]][[j]]$Palet <- bpal
            }
        }
    } else {}
    thetimer$Colour.row.factor <- (proc.time() - start)[3]
    start <- proc.time()
    if(legend){
        
        # LEGENDS
        
        cat('Creating legends...\n')
        flush.console()
        version <- installed.packages()[which(names(installed.packages()[,'Version']) =='ggplot2'), 'Version']
        
        ## Create grobs
        alegend <- list()
        
        ### Columns
        clegend <- list()
        if(!missing(colfactor) & length(colfacname)!=0){
            for(i in 1:length(colfac)){
                if(!missing(suppress.leg)){
                    if(colfacname[i] %in% suppress.leg){
                        next
                    } else {}
                } else {}
                if(colfac[[i]][[1]]$type[1] == 'seq'){
                    dat <- data.frame(x=1:9, y=1, fill=1:9)
                    l <- qplot(x, y, data=dat, geom='tile', fill=factor(fill))
                    l <- l + scale_fill_brewer(colfacname[i], breaks=1:9, labels=c(sprintf('%.4g', min(do.call('rbind', colfac[[i]])$value)), rep('', 7), sprintf('%.4g', max(do.call('rbind', colfac[[i]])$value))), palette=colfac[[i]][[1]]$Palet[1], guide=guide_legend(nrow=9))
                    l <- l + theme(legend.position='top', legend.direction='vertical')
                    if(version == '0.8.9'){
                        l <- l + opts(keep='legend_box', legend.position='left')
                    } else {
                        l <- ggplot_gtable(ggplot_build(l))
                        leg <- which(sapply(l$grobs, function(x) x$name) == "guide-box")
                        l <- l$grobs[[leg]]
                    }
                    clegend[[i]] <- l
                } else {
                    dat <- ldply(colfac[[i]], function(x) x[1,])
                    l <- qplot(x, y, data=dat, geom='tile', fill=.id)
                    l <- l + scale_fill_manual(colfacname[i], values=dat$colour, breaks=dat$.id, guide=guide_legend(nrow=9))
                    l <- l + theme(legend.position='top', legend.direction='vertical')
                    if(version == '0.8.9'){
                        l <- l + opts(keep='legend_box', legend.position='left')
                    } else {
                        l <- ggplot_gtable(ggplot_build(l))
                        leg <- which(sapply(l$grobs, function(x) x$name) == "guide-box")
                        l <- l$grobs[[leg]]
                    }
                    clegend[[i]] <- l
                }
            }
            alegend$c <- clegend
        } else {}
        
        ### Rows
        rlegend <- list()
        if(!missing(rowfactor) & length(rowfacname)!=0){
            for(i in 1:length(rowfac)){
                if(!missing(suppress.leg)){
                    if(rowfacname[i] %in% suppress.leg){
                        next
                    } else {}
                } else {}
                if(rowfac[[i]][[1]]$type[1] == 'seq'){
                    dat <- data.frame(x=1:9, y=1, fill=1:9)
                    l <- qplot(x, y, data=dat, geom='tile', fill=factor(fill))
                    l <- l + scale_fill_brewer(rowfacname[i], breaks=1:9, labels=c(sprintf('%.4g', min(do.call('rbind', rowfac[[i]])$value)), rep('', 7), sprintf('%.4g', max(do.call('rbind', rowfac[[i]])$value))), palette=rowfac[[i]][[1]]$Palet[1], guide=guide_legend(nrow=9))
                    l <- l + theme(legend.position='top', legend.direction='vertical')
                    if(version == '0.8.9'){
                        l <- l + opts(keep='legend_box', legend.position='left')
                    } else {
                        l <- ggplot_gtable(ggplot_build(l))
                        leg <- which(sapply(l$grobs, function(x) x$name) == "guide-box")
                        l <- l$grobs[[leg]]
                    }
                    rlegend[[i]] <- l
                } else {
                    dat <- ldply(rowfac[[i]], function(x) x[1,])
                    l <- qplot(x, y, data=dat, geom='tile', fill=.id)
                    l <- l + scale_fill_manual(rowfacname[i], values=dat$colour, breaks=dat$.id, guide=guide_legend(nrow=9))
                    l <- l + theme(legend.position='top', legend.direction='vertical')
                    if(version == '0.8.9'){
                        l <- l + opts(keep='legend_box', legend.position='left')
                    } else {
                        l <- ggplot_gtable(ggplot_build(l))
                        leg <- which(sapply(l$grobs, function(x) x$name) == "guide-box")
                        l <- l$grobs[[leg]]
                    }
                    rlegend[[i]] <- l
                }
            }
            alegend$r <- rlegend
        } else {}
        alegend <- do.call('c', alegend)
        
        ## Create legend plot
        cat('Plotting legends...\n')
        flush.console()
        
        pdf(file='heatLegends.pdf', height=12, width=9)
        for(i in 1:nPages){
            if(i>1) grid.newpage() else {}
            sub <- ceiling(length(alegend)/nPages)
            sub <- 1:sub + (i-1)*sub
            sub <- alegend[sub]
            do.call('grid.arrange', sub)
        }
        dev.off()
        cat(paste('heatLegend.pdf written to ', getwd(), '\n', sep=''))
    } else {
        
        # Transform intensity data
        
        cat('Transforming data...\n')
        flush.console()
        
        if(!missing(rOrder)){
            data <- data[rOrdering, ]
        } else {
            data <- data[match(rowlab, rownames(data)), ]
        }
        if(!missing(cOrder)){
            data <- data[, cOrdering]
        } else {
            data <- data[, match(collab, colnames(data))]
        }
        rownames(data) <- 1:nrow(data)
        colnames(data) <- 1:ncol(data)
        data <- data.frame(rows=row.names(data), data)
        data <- melt(data, id=1, variable.name='cols')
        data$cols <- sub('X', '', data$cols)
        data$cols <- as.numeric(data$cols)
        data$rows <- as.numeric(as.character(data$rows))
        data$gridrow <- 2+length(colfacname)
        data$gridcol <- 1
        thetimer$Transform.data <- (proc.time() - start)[3]
        start <- proc.time()
        
        # Transform dendrogram data
        
        densCoord <- data.frame(x=NA, y=NA, xend=NA, yend=NA)
        
        ## Columns
        if(missing(cOrder)){
            colsseg <- cols$segment
            colsseg$gridrow <- 1
            colsseg$gridcol <- 1
            colsseg[,c('y', 'yend')] <- colsseg[,c('y', 'yend')] + length(rowlab) + 100
            densCoord$yend <- max(colsseg[, c('y', 'yend')])
            if(!missing(cGroup)){
                colstri <- cols$triangle
                colstri$gridrow <- 1
                colstri$gridcol <- 1
                colsgtext <- colstri
                colstri$y <- colstri$y + length(rowlab) + 100
                colsgtext <- colsgtext[which(colsgtext$y != 0),]
                colsgtext$y <- min(colsgtext$y)/4 + length(rowlab) + 100
                densCoord$y <- min(colstri$y)
            } else {
                densCoord$y <- min(colsseg[, c('y', 'yend')])
            }
        } else {
            if(cOrdertype == 'qual'){
                dencolheight <- max(cols$segment$y)
                colstri <- ddply(cOrderlev, .(order), function(x) data.frame(x=c(x$min, x$max, mean(c(x$min, x$max))), y=c(0,0,dencolheight), group=x$factor))
                colstri$gridrow <- 1
                colstri$gridcol <- 1
                colsgtext <- colstri
                colstri$y <- colstri$y + length(rowlab) + 100
                colsgtext <- colsgtext[which(colsgtext$y != 0),]
                colsgtext$y <- min(colsgtext$y)/6 + length(rowlab) + 100
                fromto <- data.frame(from=1:nrow(colsgtext), to=colsgtext$x)
                for(i in nrow(fromto):1){
                    cols <- denRescale(cols, fromto[i, 'from'], fromto[i, 'to'])
                }
                colsseg <- cols$segment
                colsseg$gridrow <- 1
                colsseg$gridcol <- 1
                colsseg[,c('y', 'yend')] <- colsseg[,c('y', 'yend')] + length(rowlab) + 100 + dencolheight
                densCoord$yend <- max(colsseg[, c('y', 'yend')])
                densCoord$y <- min(colstri$y)
            } else {
                carrow <- data.frame(x=(1 + length(collab)/10), xend=(length(collab)*0.9), y=length(rowlab) + 100, yend=length(rowlab) + 100, gridrow=1, gridcol=1)
                climit <- rbind(carrow, carrow)
                climit[1, c('y', 'yend')] <- climit[1, c('y', 'yend')] - 1
                climit[2, c('y', 'yend')] <- climit[2, c('y', 'yend')] + 2
                carrow.title <- data.frame(title=paste('Sorted by ', cOrder[1], sep=''), x=mean(c(1, length(collab))), y=length(rowlab) + 100, gridrow=1, gridcol=1)
                densCoord$y <- climit$y[1]
                densCoord$yend <- climit$y[2]
            }
        }
        
        ## Rows
        if(missing(rOrder)){
            rowsseg <- rows$segment
            rowsseg$gridrow <- 2+length(colfacname)
            rowsseg$gridcol <- 2+length(rowfacname)
            rowsseg[,c('y', 'yend')] <- rowsseg[,c('y', 'yend')] + length(collab) + 100
            densCoord$xend <- max(rowsseg[, c('y', 'yend')])
            if(!missing(rGroup)){
                rowstri <- rows$triangle
                rowstri$gridrow <- 2+length(colfacname)
                rowstri$gridcol <- 2+length(rowfacname)
                rowsgtext <- rowstri
                rowstri$y <- rowstri$y + length(collab) + 100
                rowsgtext <- rowsgtext[which(rowsgtext$y != 0),]
                rowsgtext$y <- min(rowsgtext$y)/2 + length(collab) + 100
                densCoord$x <- min(rowstri$y)
            } else {
                densCoord$x <- min(rowsseg[, c('y', 'yend')])
            }
        } else {
            if(rOrdertype == 'qual'){
                denrowheight <- max(rows$segment$y)
                rowstri <- ddply(rOrderlev, .(order), function(x) data.frame(x=c(x$min, x$max, mean(c(x$min, x$max))), y=c(0,0,denrowheight), group=x$factor))
                rowstri$gridrow <- 2+length(colfacname)
                rowstri$gridcol <- 2+length(rowfacname)
                rowsgtext <- rowstri
                rowstri$y <- rowstri$y + length(collab) + 100
                rowsgtext <- rowsgtext[which(rowsgtext$y != 0),]
                rowsgtext$y <- min(rowsgtext$y)/50 + length(collab) + 100
                fromto <- data.frame(from=1:nrow(rowsgtext), to=rowsgtext$x)
                for(i in nrow(fromto):1){
                    rows <- denRescale(rows, fromto[i, 'from'], fromto[i, 'to'])
                }
                rowsseg <- rows$segment
                rowsseg$gridrow <- 2+length(colfacname)
                rowsseg$gridcol <- 2+length(rowfacname)
                rowsseg[,c('y', 'yend')] <- rowsseg[,c('y', 'yend')] + length(collab) + 100 + denrowheight
                densCoord$xend <- max(rowsseg[, c('y', 'yend')])
                densCoord$x <- min(rowstri$y)
            } else {
                rarrow <- data.frame(x=(1 + length(rowlab)/10), xend=(length(rowlab)*0.9), y=length(collab) + 100, yend=length(collab) + 100, gridrow=2+length(colfacname), gridcol=2+length(rowfacname))
                rlimit <- rbind(rarrow, rarrow)
                rlimit[1, c('y', 'yend')] <- rlimit[1, c('y', 'yend')] - 1
                rlimit[2, c('y', 'yend')] <- rlimit[2, c('y', 'yend')] + 2
                rarrow.title <- data.frame(title=paste('Sorted by ', rOrder[1], sep=''), x=mean(c(1, length(rowlab))), y=length(collab) + 100, gridrow=2+length(colfacname), gridcol=2+length(rowfacname))
                densCoord$x <- rlimit$y[1]
                densCoord$xend <- rlimit$y[2]
            }
        }
        thetimer$Transform.dendrogram <- (proc.time() - start)[3]
        start <- proc.time()
        
        # Create density distribution data
        
        dens <- density(data$value)
        dens <- data.frame(x=dens$x, y=dens$y, gridrow=1, gridcol=2+length(rowfacname))
        dens <- dens[which(dens$x >= range(data$value)[1] & dens$x <= range(data$value)[2]),]
        minmaxr <- as.numeric(densCoord[, c('x', 'xend')])
        minmaxr <- minmaxr + c(diff(minmaxr)/4, -diff(minmaxr)/4)
        dens$x <- dens$x*diff(minmaxr)/max(dens$x) + minmaxr[1]
        dens$y <- dens$y + (0 - min(dens$y))
        minmaxc <- as.numeric(densCoord[, c('y', 'yend')])
        minmaxc <- minmaxc + c(diff(minmaxc)/4, -diff(minmaxc)/4)
        dens$y <- dens$y*diff(minmaxc)/max(dens$y) + minmaxc[1]
        
        # Create background for density plot
        
        x0 <- seq(-0.5, 9.5)*diff(minmaxr)/10 + minmaxr[1]
        x1 <- seq(0.5, 10.5)*diff(minmaxr)/10 + minmaxr[1]
        densback <- data.frame(x0=x0, x1=x1, y0=minmaxc[1]-diff(minmaxc)*0.05, y1=minmaxc[2]+diff(minmaxc)*0.05, value=seq(0, 100, by=10), gridrow=1, gridcol=2+length(rowfacname))
        thetimer$Create.density <- (proc.time() - start)[3]
        start <- proc.time()
        
        # Restore dimnames
        
        collab <- rawDimNames$col[as.numeric(collab)]
        rowlab <- rawDimNames$row[as.numeric(rowlab)]
        
        # Create plots
        
        cat('Creating plots...\n')
        flush.console()
        
        ## Heatmap
        data$gridcol <- factor(data$gridcol, levels=1:(2+length(rowfacname)))
        data$gridrow <- factor(data$gridrow, levels=1:(2+length(colfacname)))
        p <- ggplot(data=data) + theme_bw()
        p <- p + geom_raster(aes(x=cols, y=rows, fill=value))
        
        ## Factors
        
        ### For columns
        if(!missing(colfactor) & length(colfacname)!=0){
            for(i in 1:length(colfac)){
                for(j in 1:length(colfac[[i]])){
                    colfac[[i]][[j]]$gridcol <- factor(colfac[[i]][[j]]$gridcol, levels=1:(2+length(rowfacname)))
                    colfac[[i]][[j]]$gridrow <- factor(colfac[[i]][[j]]$gridrow, levels=1:(2+length(colfacname)))
                    p <- p + geom_tile(data=colfac[[i]][[j]], aes(x=x, y=y, fill=NULL, width=1), fill=colfac[[i]][[j]]$colour[1])
                }
            }
        } else {}
        
        ### for Rows
        if(!missing(rowfactor) & length(rowfacname)!=0){
            for(i in 1:length(rowfac)){
                for(j in 1:length(rowfac[[i]])){
                    rowfac[[i]][[j]]$gridcol <- factor(rowfac[[i]][[j]]$gridcol, levels=1:(2+length(rowfacname)))
                    rowfac[[i]][[j]]$gridrow <- factor(rowfac[[i]][[j]]$gridrow, levels=1:(2+length(colfacname)))
                    p <- p + geom_tile(data=rowfac[[i]][[j]], aes(x=x, y=y, fill=NULL, height=1), fill=rowfac[[i]][[j]]$colour[1])
                }
            }
        } else {}
        thetimer$Build.factors <- (proc.time() - start)[3]
        start <- proc.time()
        
        ## Dendrograms
        
        ### Columns
        if(cOrdertype != 'seq'){
            colsseg$gridcol <- factor(colsseg$gridcol, levels=1:(2+length(rowfacname)))
            colsseg$gridrow <- factor(colsseg$gridrow, levels=1:(2+length(colfacname)))
            p <- p + geom_segment(data=colsseg, aes(x=x, y=y, xend=xend, yend=yend))
            if(!missing(cGroup) | !missing(cOrder)){
                colstri$gridcol <- factor(colstri$gridcol, levels=1:(2+length(rowfacname)))
                colstri$gridrow <- factor(colstri$gridrow, levels=1:(2+length(colfacname)))
                colsgtext$gridcol <- factor(colsgtext$gridcol, levels=1:(2+length(rowfacname)))
                colsgtext$gridrow <- factor(colsgtext$gridrow, levels=1:(2+length(colfacname)))
                p <- p + geom_polygon(data=colstri, aes(x=x, y=y, group=group), colour='black', fill='grey80')
                p <- p + geom_text(data=colsgtext, aes(x=x, y=y, label=group), size=3)
            } else {}	  	
        } else {
            carrow$gridcol <- factor(carrow$gridcol, levels=1:(2+length(rowfacname)))
            carrow$gridrow <- factor(carrow$gridrow, levels=1:(2+length(colfacname)))
            climit$gridcol <- factor(climit$gridcol, levels=1:(2+length(rowfacname)))
            climit$gridrow <- factor(climit$gridrow, levels=1:(2+length(colfacname)))
            carrow.title$gridcol <- factor(carrow.title$gridcol, levels=1:(2+length(rowfacname)))
            carrow.title$gridrow <- factor(carrow.title$gridrow, levels=1:(2+length(colfacname)))
            p <- p + geom_segment(data=carrow, aes(x=x, y=y, xend=xend, yend=yend), arrow=arrow(length=unit(0.5, 'cm')))
            p <- p + geom_segment(data=climit, aes(x=x, y=y, xend=xend, yend=yend), alpha=(0))
            p <- p + geom_text(data=carrow.title, aes(x=x, y=y, label=title), vjust=-0.5, size=4)
        }
        
        ### Rows
        if(rOrdertype != 'seq'){
            rowsseg$gridcol <- factor(rowsseg$gridcol, levels=1:(2+length(rowfacname)))
            rowsseg$gridrow <- factor(rowsseg$gridrow, levels=1:(2+length(colfacname)))
            p <- p + geom_segment(data=rowsseg, aes(x=y, y=x, xend=yend, yend=xend))
            if(!missing(rGroup) | !missing(rOrder)){
                rowstri$gridcol <- factor(rowstri$gridcol, levels=1:(2+length(rowfacname)))
                rowstri$gridrow <- factor(rowstri$gridrow, levels=1:(2+length(colfacname)))
                rowsgtext$gridcol <- factor(rowsgtext$gridcol, levels=1:(2+length(rowfacname)))
                rowsgtext$gridrow <- factor(rowsgtext$gridrow, levels=1:(2+length(colfacname)))
                p <- p + geom_polygon(data=rowstri, aes(x=y, y=x, group=group), colour='black', fill='grey80')
                p <- p + geom_text(data=rowsgtext, aes(x=y, y=x, label=group), size=3, hjust=0)
            } else {}
        } else {
            rarrow$gridcol <- factor(rarrow$gridcol, levels=1:(2+length(rowfacname)))
            rarrow$gridrow <- factor(rarrow$gridrow, levels=1:(2+length(colfacname)))
            rlimit$gridcol <- factor(rlimit$gridcol, levels=1:(2+length(rowfacname)))
            rlimit$gridrow <- factor(rlimit$gridrow, levels=1:(2+length(colfacname)))
            rarrow.title$gridcol <- factor(rarrow.title$gridcol, levels=1:(2+length(rowfacname)))
            rarrow.title$gridrow <- factor(rarrow.title$gridrow, levels=1:(2+length(colfacname)))
            p <- p + geom_segment(data=rarrow, aes(x=y, y=x, xend=yend, yend=xend), arrow=arrow(length=unit(0.5, 'cm')))
            p <- p + geom_segment(data=rlimit, aes(x=y, y=x, xend=yend, yend=xend), alpha=(0))
            p <- p + geom_text(data=rarrow.title, aes(x=y, y=x, label=title), vjust=-0.5, size=4, angle=-90)
        }
        thetimer$Build.dendro <- (proc.time() - start)[3]
        start <- proc.time()
        
        ## Distribution
        densback$gridcol <- factor(densback$gridcol, levels=1:(2+length(rowfacname)))
        densback$gridrow <- factor(densback$gridrow, levels=1:(2+length(colfacname)))
        dens$gridcol <- factor(dens$gridcol, levels=1:(2+length(rowfacname)))
        dens$gridrow <- factor(dens$gridrow, levels=1:(2+length(colfacname)))
        p <- p + geom_rect(data=densback, aes(xmin=x0, xmax=x1, ymin=y0, ymax=y1, fill=value))
        p <- p + geom_line(data=dens, aes(x=x, y=y), colour='white', size=I(2))
        thetimer$Build.distribution <- (proc.time() - start)[3]
        start <- proc.time()
        
        ## Faceting
        p <- p + facet_grid(gridrow~gridcol, scales='free', height=c(1-0.05*length(colfacname),rep(0.05, length(colfacname)),2), width=c(2,rep(0.05, length(rowfacname)),1-0.05*length(rowfacname)))
        
        ## Set scales
        if('columns' %in% hideticks){
            if(!missing(rowfactor) & length(rowfacname)!=0){
                p <- p + scale_x_continuous('', breaks=-1000*1:length(rowfacname), labels=rowfacname, expand=c(0,0))
            } else {
                p <- p + scale_x_continuous('', breaks=-1, labels='', expand=c(0,0))
            }
        } else {
            if(!missing(rowfactor) & length(rowfacname)!=0){
                p <- p + scale_x_continuous('', breaks=c(-1000*1:length(rowfacname), 1:length(collab)), labels=c(rowfacname, collab), expand=c(0,0))
            } else {
                p <- p + scale_x_continuous('', breaks=1:length(collab), labels=collab, expand=c(0,0))
            }
        }
        if('rows' %in% hideticks){
            if(!missing(colfactor) & length(colfacname)!=0){
                p <- p + scale_y_continuous('', breaks=-1000*1:length(colfacname), labels=colfacname, expand=c(0,0))
            } else {
                p <- p + scale_y_continuous('', breaks=-1, labels='', expand=c(0,0))
            }
        } else {
            if(!missing(colfactor) & length(colfacname)!=0){
                p <- p + scale_y_continuous('', breaks=c(-1000*1:length(colfacname), 1:length(rowlab)), labels=c(colfacname, rowlab), expand=c(0,0))
            } else {
                p <- p + scale_y_continuous('', breaks=1:length(rowlab), labels=rowlab, expand=c(0,0))
            }
        }
        p <- p + scale_fill_gradientn(dataname, colours=c('blue', 'black', 'red'), space='Lab', guide=guide_colourbar())
        
        ## Set options
        p <- p + theme(axis.text.x=element_text(angle=-90, hjust=0, size=8, vjust=0.2), axis.text.y=element_text(hjust=1, size=8))
        p <- p + theme(panel.border=element_blank()) + theme(strip.text.y=element_blank(), strip.text.x=element_blank(), strip.background=element_blank())
        p <- p + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
        p <- p + ggtitle(title)
        thetimer$Setup.plot <- (proc.time() - start)[3]
        start <- proc.time()
        
        # Create output
        
        cat('Plotting...\n')
        flush.console()
        
        ans <- list()
        if(!missing(cGroup)){
            names(cGroup) <- collab
            ans$Columns <- cGroup
        } else {}
        if(!missing(rGroup)){
            names(rGroup) <- rowlab
            ans$Rows <- rGroup
        } else {}
        if(outputPlot | length(ans) == 0){
            p
        } else {
            print(p)
            invisible(ans)
        }
    }
}
### Grid layout functions
vplayout <- function (...) {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(...)))
}
subplot <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)