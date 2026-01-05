gc.sample.stats <- function(file, col_types = "c--dd----d----",
    buffer = 335544032, parallel = 2L, verbose = TRUE) {
    con <- gzfile(file, "rb")
    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, col_types) {
        x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, col_names = FALSE,
            skip = 0, n_max = Inf, progress = FALSE)
        u_chr <- unique(x[, 1])
        n_chr <- table(x[, 1])
        gc1 <- lapply(split(x[, 2], x[, 4]), table)
        gc2 <- lapply(split(x[, 3], x[, 4]), table)
        if (verbose){
            message(".", appendLF = FALSE)
        }
        list(unique = u_chr, lines = n_chr, gc_nor = gc1, gc_tum = gc2)
    }
    if (verbose){
        message("Collecting GC information ", appendLF = FALSE)
    }
    res <- chunk.apply(input = con, FUN = parse_chunck, col_types = col_types,
        CH.MAX.SIZE = buffer)
    close(con)
    if (verbose){
        message(" done\n")
    }
    unfold_gc(res, stats = TRUE)
}
unfold_gc <- function(x, stats = TRUE) {
    gc_norm <- get_gc(x[, "gc_nor"])
    gc_tum <- get_gc(x[, "gc_tum"])
    if (stats) {
        ord_chrom <- unique(Reduce("c", Reduce("c", x[, "unique"])))
        stats_chrom <- Reduce("c", x[, "lines"])
        stats_chrom <- sapply(splash_table(x[, "lines"]), sum)
        stats_chrom <- stats_chrom[ord_chrom]
        stats_start <- cumsum(c(1, stats_chrom[-length(stats_chrom)]))
        stats_end   <- stats_start + stats_chrom - 1
        stats_chrom <- data.frame(chr = ord_chrom, n_lines = stats_chrom,
            start = stats_start, end = stats_end)
        list(file.metrics = stats_chrom, normal = gc_norm, tumor = gc_tum)
    } else {
        list(normal = gc_norm, tumor = gc_tum)
    }
}
splash_table <- function(lis_obj){
    lis_obj <- Reduce("c", lis_obj)
    split(lis_obj, names(lis_obj))
}
get_gc <- function(gc_col) {
    sort_char <- function(x) {
        as.character(sort(as.numeric(x)))
    }
    all_depths <- splash_table(gc_col)
    all_depths <- lapply(all_depths, FUN = function(x) {
        sapply(splash_table(x), sum)
    })
    names_gc <- sort_char(names(all_depths))
    all_depths <- all_depths[names_gc]
    names_depths <- sort_char(unique(Reduce("c", lapply(all_depths, names))))
    n <- do.call(rbind, lapply(all_depths, FUN = function(x, names_depths) {
            res <- x[names_depths]
            names(res) <- names_depths
            res
        },
        names_depths = names_depths))
    n[is.na(n)] <- 0
    list(gc = as.numeric(names_gc), depth = as.numeric(names_depths), n = n)
}
median_gc <- function(gc_list) {
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.median(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}
mean_gc <- function(gc_list) {
    apply(gc_list$n, 1, FUN = function(x, w) {
            weighted.mean(x = w, w = x, na.rm = TRUE)
        },
        w = gc_list$depth)
}
depths_gc <- function(depth_n, depth_t, gc) {
    gc_nor <- lapply(split(depth_n, gc), table)
    gc_tum <- lapply(split(depth_t, gc), table)
    list(gc_nor = gc_nor, gc_tum = gc_tum)
}
read.seqz <- function(file, n_lines = NULL, col_types = "ciciidddcddccc",
    chr_name = NULL, buffer = 33554432, parallel = 1,
    col_names = c("chromosome", "position", "base.ref", "depth.normal",
                  "depth.tumor", "depth.ratio", "Af", "Bf", "zygosity.normal",
                  "GC.percent", "good.reads", "AB.normal", "AB.tumor",
                  "tumor.strand"), ...) {

    if (is.null(n_lines)) {
        skip <- 1
        n_max <- Inf
    } else {
        n_lines <- round(sort(n_lines), 0)
        skip <- n_lines[1]
        n_max <- n_lines[2] - skip + 1
    }
    if (!is.null(chr_name)) {
        chr_name <- as.character(chr_name)
        tbi <- file.exists(paste(file, "tbi", sep = "."))
        if(FALSE){
            read.seqz.tbi(file, split_chr_coord(chr_name),
                col_names, col_types)
        } else {
            read.seqz.chr(file, chr_name = chr_name, col_types = col_types,
                col_names = col_names, buffer = buffer, parallel = parallel)
        }
    } else {
        read_tsv(file = file, col_types = col_types, skip = skip,
            n_max = n_max, col_names = col_names, progress = FALSE, ...)
    }
}
read.seqz.chr <- function(file, chr_name, col_names,
    col_types, buffer, parallel) {
    con <- gzfile(file, "rb")
    suppressWarnings(skip_line <- readLines(con, n = 1))
    remove(skip_line)
    parse_chunck <- function(x, chr_name, col_names, col_types) {
        x <- read_tsv(file = paste(mstrsplit(x), collapse = "\n"),
            col_types = col_types, skip = 0, n_max = Inf,
            col_names = col_names, progress = FALSE)
        x[x$chromosome == chr_name, ]
    }
    res <- chunk.apply(input = con, FUN = parse_chunck, chr_name = chr_name,
        col_names = col_names, col_types = col_types, CH.MAX.SIZE = buffer)
    close(con)
    res
}
read.seqz.tbi <- function(file, chr_name, col_names, col_types) {
    res <- tabix.read(file, chr_name)
    res <- read_tsv(file = paste(mstrsplit(res), collapse = "\n"),
        col_types = col_types, skip = 0, n_max = Inf,
        col_names = col_names, progress = FALSE)
}
sequenza.extract_custom <- function(file, window = 1e4, overlap = 1,
    gamma = 80, kmin = 10, gamma.pcf = 140, kmin.pcf = 40,
    mufreq.treshold = 0.10, min.reads = 40, min.reads.normal = 10,
    min.reads.baf = 1, max.mut.types = 1, min.type.freq = 0.9,
    min.fw.freq = 0, verbose = TRUE, chromosome.list = NULL,
    breaks = NULL, breaks.method = "het", assembly = "hg19",
    weighted.mean = TRUE, normalization.method = "mean",
    ignore.normal = FALSE, parallel = 1, gc.stats = NULL,
    segments.samples = FALSE){
    if (is.null(gc.stats)) {
        gc.stats <- gc.sample.stats(file, verbose = verbose,
            parallel = parallel)
    }
    if (normalization.method == "mean") {
        gc.normal.vect <- mean_gc(gc.stats$normal)
        gc.tumor.vect  <- mean_gc(gc.stats$tumor)
        avg_tum_depth <- weighted.mean(x = gc.stats$tumor$depth,
            w = colSums(gc.stats$tumor$n))
        avg_nor_depth <- weighted.mean(x = gc.stats$normal$depth,
            w = colSums(gc.stats$normal$n))
    } else {
        gc.normal.vect <- median_gc(gc.stats$normal)
        gc.tumor.vect  <- median_gc(gc.stats$tumor)
        avg_tum_depth <- weighted.median(x = gc.stats$tumor$depth,
            w = colSums(gc.stats$tumor$n))
        avg_nor_depth <- weighted.median(x = gc.stats$normal$depth,
            w = colSums(gc.stats$normal$n))
    }

    segments_samples.list <- list()
    norm.gc.list <- list()

    seqz_hetero_list<-list();
    seqz.list<-list();

    if (is.null(dim(breaks))) {
        breaks <- NULL
    }
    chr.vect <- as.character(gc.stats$file.metrics$chr)
    if (is.null(chromosome.list)) {
        chromosome.list <- chr.vect
    } else {
        chromosome.list <- chromosome.list[chromosome.list %in% chr.vect]
    }




    for (chr in chromosome.list) {
        if (verbose) {
            message("Processing ", chr, ":", appendLF = TRUE)
        }
        tbi <- file.exists(paste0(file, ".tbi") )
        if (tbi) {
            seqz.data   <- read.seqz(file, chr_name = chr)
        } else {
            file.lines <- gc.stats$file.metrics[which(chr.vect == chr), ]
            seqz.data   <- read.seqz(file, n_lines = c(file.lines$start,
                file.lines$end))
        }
        norm_tumor_depth <- seqz.data$depth.tumor /
            gc.tumor.vect[as.character(seqz.data$GC.percent)]
        norm_normal_depth <- seqz.data$depth.normal /
            gc.normal.vect[as.character(seqz.data$GC.percent)]


        norm.gc.stats <- depths_gc(
            depth_n = round(norm_normal_depth * avg_nor_depth, 0),
            depth_t = round(norm_tumor_depth * avg_tum_depth, 0),
            gc = seqz.data$GC.percent)
        if (ignore.normal) {
            seqz.data$adjusted.ratio <- round(norm_tumor_depth, 3)
        } else {
            seqz.data$adjusted.ratio <- round(
                norm_tumor_depth / norm_normal_depth, 3)
        }



        seqz.hom <- seqz.data$zygosity.normal == "hom"
        seqz.het <- seqz.data[!seqz.hom, ]
        het.filt <- seqz.het$good.reads >= min.reads.baf
        seqz.het <- seqz.het[het.filt, ]
        het_ok <- nrow(seqz.het) > 0

        if (verbose) {
            message('   ', nrow(seqz.het), ' heterozygous positions.', appendLF = TRUE)
            message('   ', sum(seqz.hom), ' homozygous positions.', appendLF = TRUE)
        }

        norm.gc.list[[chr]] <- norm.gc.stats


        seqz.list[[chr]]<-seqz.data
        seqz_hetero_list[[chr]]<-seqz.het;
    }

    gc_norm <- unfold_gc(do.call(rbind, norm.gc.list), stats = FALSE)

    avg_tum_ndepth <- weighted.mean(x = gc_norm$tumor$depth,
                                    w = colSums(gc_norm$tumor$n))
    avg_nor_ndepth <- weighted.mean(x = gc_norm$normal$depth,
                                    w = colSums(gc_norm$normal$n))
    avg_depth_ratio <- (avg_tum_ndepth / avg_tum_depth) /
        (avg_nor_ndepth / avg_nor_depth)

    print(paste0("avg_depth_ratio=",avg_depth_ratio))
    list(
        b=seqz.list,
        d=seqz_hetero_list,
        avg_depth_ratio=avg_depth_ratio)
}
segment.breaks <- function(seqz.tab, breaks, min.reads.baf = 1,
    weighted.mean = FALSE) {
    if (weighted.mean){
        w.r <- sqrt(seqz.tab$depth.normal)
        rw <- seqz.tab$adjusted.ratio * w.r
        w.b <- sqrt(seqz.tab$good.reads)
        bw <- seqz.tab$Bf * w.b
        seqz.tab <- cbind(seqz.tab[, c("chromosome", "position",
            "zygosity.normal", "good.reads", "Af", "Bf")],
            rw = rw, w.r = w.r, bw = bw, w.b = w.b)
    }
    chr.order <- unique(seqz.tab$chromosome)
    seqz.tab <- split(seqz.tab, f = seqz.tab$chromosome)
    segments <- list()
    for (i in 1:length(seqz.tab)) {
        seqz.b.i <- seqz.tab[[i]][seqz.tab[[i]]$zygosity.normal == "het", ]
        seqz.b.i <- seqz.b.i[seqz.b.i$good.reads >= min.reads.baf, ]
        breaks.i <- breaks[breaks$chrom == names(seqz.tab)[i], ]
        nb <- nrow(breaks.i)
        breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,
            c("start.pos", "end.pos")], f = 1:nb))
        unique.breaks <- function(b, offset = 1) {
            while(any(diff(b) == 0)) {
                b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + offset
            }
            b
        }
        breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), offset = 1)
        fact.r.i <- cut(seqz.tab[[i]]$position, breaks.vect)
        fact.b.i <- cut(seqz.b.i$position, breaks.vect)
        seg.i.s.r <- sapply(X = split(seqz.tab[[i]]$chromosome,
            f = fact.r.i), FUN = length)
        seg.i.s.b <- sapply(X = split(seqz.b.i$chromosome,
            f = fact.b.i), FUN = length)
        if (weighted.mean) {
            seg.i.rw    <- sapply(X = split(seqz.tab[[i]]$rw, f = fact.r.i),
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.w.r   <- sapply(X = split(seqz.tab[[i]]$w.r, f = fact.r.i),
                FUN = function(a) sum(a, na.rm = TRUE))
            seg.i.r.sd  <- sapply(X = split(seqz.tab[[i]]$rw /
                seqz.tab[[i]]$w.r, f = fact.r.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd  <- sapply(X = split(seqz.b.i$bw /
                seqz.b.i$w.b, f = fact.b.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            A.split <- split(seqz.b.i$Af, f = fact.b.i)
            B.split <- split(seqz.b.i$Bf, f = fact.b.i)
            d.split <- split(seqz.b.i$good.reads, f = fact.b.i)
            window.quantiles <- mapply(b_allele_freq, Af = A.split,
                Bf = B.split, good.reads = d.split, conf = 0.95)
            segments.i <- data.frame(chromosome  = names(seqz.tab)[i],
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                end.pos = as.numeric(breaks.vect[-1]),
                Bf = window.quantiles[2, ], N.BAF = seg.i.s.b,
                sd.BAF = seg.i.b.sd, depth.ratio = seg.i.rw / seg.i.w.r,
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
                stringsAsFactors = FALSE)
        } else {
            seg.i.r    <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio,
                f = fact.r.i), FUN = function(a) mean(a, na.rm = TRUE))
            A.split <- split(seqz.b.i$Af, f = fact.b.i)
            B.split <- split(seqz.b.i$Bf, f = fact.b.i)
            d.split <- split(seqz.b.i$good.reads, f = fact.b.i)



            X = split(seqz.tab[[i]]$Bf, f = fact.r.i)
            Z = split(seqz.tab[[i]]$zygosity.normal, f = fact.r.i);

            seg.i.b <- mapply(function(x, z){mean(x[z=="het"], na.rm=TRUE) }, x = X,
                                       z = Z)




            seg.i.r.sd <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio,
                f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
            seg.i.b.sd <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i),
                FUN = function(a) sd(a, na.rm = TRUE))
            segments.i <- data.frame(chromosome  = names(seqz.tab)[i],
                start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                end.pos = as.numeric(breaks.vect[-1]),
                Bf = seg.i.b, N.BAF = seg.i.s.b,
                sd.BAF = seg.i.b.sd, depth.ratio = seg.i.r,
                N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd,
                stringsAsFactors = FALSE)
        }
        segments[[i]] <- segments.i[seq(from = 1,
            to = nrow(segments.i), by = 2),]
    }
    segments <- do.call(rbind, segments[as.factor(chr.order)])
    row.names(segments) <- 1:nrow(segments)
    len.seg <- (segments$end.pos - segments$start.pos) / 1e6
    segments[(segments$N.ratio / len.seg) >= 2, ]
}
dt2 <- function(x, df, ncp, log = FALSE, mean, sd) {
  x2 <- (x - mean) / sd
  dt(x2, df = df, ncp = ncp, log = log)
}
theoretical.depth.ratio2 <- function(cellularity, B,  ploidy, CNt,CNn = 2,
                                    normal.ploidy = 2, avg.depth.ratio = 1) {


  cellu_copy_term <- (1 - cellularity) + (CNt / CNn * cellularity)
  ploidy_cellu_term <- (ploidy / normal.ploidy * cellularity) +
    1 - cellularity

  re<-avg.depth.ratio * cellu_copy_term / ploidy_cellu_term

  re

}
theoretical.baf2 <- function(CNt, B, cellularity, CNn = 2, sd.Bf=sd.Bf) {
  baf <- ( (B * cellularity) + ( 1 - cellularity) ) /
    ( (CNt * cellularity) + CNn * ( 1 - cellularity) )
  baf[CNn <= 1] <- NA
  baf
}
baf.types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
  cn_ratio_vect <- seq(from = CNt.min / CNn, to = CNt.max / CNn,
                       by = 1 / CNn)
  CNt <- cn_ratio_vect * CNn
  if (CNn < 2) {
    b_comb <- lapply(CNt, FUN = function(x) 0)
  } else {
    b_comb <- lapply(CNt, FUN = function(x) {
      seq(from = 0, to = trunc(x / 2))
    })
  }
  times_b <- sapply(b_comb, length)
  CNt <- rep(CNt, times = times_b)
  B <- unlist(b_comb)
  data.frame(CNn = CNn, CNt = CNt, B = B)
}
baf.model.points_h <- function(cellularity, ploidy,
                               avg.depth.ratio,
                               sd.Bf=sd.Bf, sd.ratio=sd.ratio,
                               t_CNt.max=t_CNt.max) {

  baf_types <- baf.types.matrix(CNt.min =1,
                                CNt.max = t_CNt.max, CNn = 2);

  depth_ratio <- theoretical.depth.ratio2(cellularity = cellularity,
                                          B = baf_types[, "B"],
                                         ploidy = ploidy,
                                         CNn = baf_types[, "CNn"],
                                         CNt = baf_types[, "CNt"],
                                         avg.depth.ratio = avg.depth.ratio);

  trunca_ratio_t<-log2( pt(((2.5)-(depth_ratio))/sd.ratio,df=5 ) - pt((0-depth_ratio)/sd.ratio,df=5 ) );


  baf <- theoretical.baf2(cellularity = cellularity,
                          CNn = baf_types[, "CNn"],
                        CNt = baf_types[, "CNt"],
                        B = baf_types[, "B"],
                        sd.Bf=sd.Bf);


  trunca_baf_t<-log2(pt((0.5-baf)/sd.Bf,df=5 ) - pt((0-baf)/sd.Bf,df=5 ) );



  model_baf<-data.frame(BAF = baf,
             depth.ratio = depth_ratio,
             trunca_ratio_t=trunca_ratio_t,
             trunca_baf_t=trunca_baf_t)


  model.pts <- cbind(CNt =  baf_types$CNt, A = baf_types$CNt - baf_types$B,
                     B = baf_types$B, model_baf);
  model.pts
}
