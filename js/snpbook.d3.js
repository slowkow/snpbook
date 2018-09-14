function plotQuery(query) {
  if (query.indexOf(':') > -1) {
    update(parseRegion(query))
  } else {
    // if (query.indexOf('rs') == 0) {
    myVariant(query)
  }
}

// TODO: Use ES6 Promises to untangle this chain of callbacks.
function update(obj) {

  console.log("[update]")
  console.log(obj)

  //var limit = -1
  var limit = 100

  var target_id = '#main-plot'

  // var data_url = 'chr6.1kg.phase3.v5a.vcf.gz.txt'
  var data_url = dataURL(obj.chrom, obj.start, obj.end, limit)

  d3.select("#loading-text").text("Loading...")

  console.log("downloading genotypes...")
  d3.xhr(data_url, "text/plain", function (err, response) {

    console.log("calculating and plotting...")
    updateLD(response, obj)
    plotLD(target_id)

    updateGenes(obj, function() { return plotGenes(target_id) })

    d3.select("#loading-text").text("")
  })

}

function myVariant(query) {
  // TODO:
  // 1. Use the iobio service.
  // 2. Handle other VCF files like ExAC.
  //
  // var data_url = "http://tabix.iobio.io/?cmd=-h%20%27http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/vcf.b37/chr2.1kg.phase3.v5a.vcf.gz%27%202:4000000-4050000"

  var data_url = 'https://myvariant.info/v1/query?q=' +
    query + '&fields=dbsnp.hg19,dbsnp.rsid,dbsnp.chrom,dbsnp.alt,dbsnp.ref'

  d3.json(data_url, function(err, json) {
    if (json.hits.length) {
      // Take the first (best) hit.
      var h = json.hits[0]
      var chrom = h.dbsnp.chrom,
          start = h.dbsnp.hg19.start,
          end = h.dbsnp.hg19.end,
          rsid = h.dbsnp.rsid
      // Update the region input.
      var el = document.getElementById('region')
      var cf = d3.format(",")
      el.value = chrom + ':' + cf(start) + '-' + cf(end)
      // Update the plot.
      update({
        chrom: chrom,
        start: start - 50000,
        end: end + 50000,
        pos: start,
        rsid: rsid
      })
    } else {
      d3.select("#loading-text")
        .text("[myvariant.info] No matches for " + query)
    }
  })
}

// Fills a global called 'data.markers'.
function updateLD(response, obj) {
    // Exclude the header lines.
    var rows = d3.tsv.parseRows(
      response.responseText.split('\n').filter(
        function(x) { return x[0] != "#" }
      ).join('\n')
    )

    // Discard multi-allelic markers.
    rows = rows.filter(function(x) { return x[4].indexOf(',') == -1 })

    // Convert data lines to a genotype matrix.
    // var geno = rows.map(get_geno)
    var geno = rows.map(genotypes)

    // Compute the minor allele frequencies.
    var maf = geno.map(function(y) {
        var retval = d3.sum(y) / y.length
        return Math.min(retval, 1 - retval)
    })

    // Choose a reference marker if one is given.
    var idx = Math.floor(geno.length / 2)
    if (obj.rsid) {
      var i = rows.map(function(x) {return x[2]}).indexOf(obj.rsid)
      if (i > -1) {
        idx = i
        console.log("Focusing on " + obj.rsid + " at index " + idx)
      }
    }

    var ld = geno.map(function(y) {
      return compute_ld(geno[idx], y)
    })

    var dist = rows.map(function(x) {
      return Math.abs(rows[idx][1] - x[1])
    })

    // // Compute D' (D prime).
    // var dp = geno.map(function(y) { 
    //   return dprime(geno[idx], y)
    // })

    // // Compute the squared Pearson correlations between centered markers.
    // var r2 = geno.map(function(y) { 
    //   var r = pearson(geno[idx], y)
    //   return r * r
    // })

    // Add the markers in this genomic region.
    data.markers = new Array(rows.length)

    for (var i = 0; i < data.markers.length; i++) {
      data.markers[i] = {
        'dist': dist[i],
        'r2': ld[i].r2,
        'dp': ld[i].dp,
        'maf': maf[i],
        'chrom': rows[i][0],
        'pos': +rows[i][1],
        'rsid': rows[i][2],
        'ref': rows[i][3],
        'alt': rows[i][4],
        'chosen': i == idx ? rows[i][2] : '',
        'geno': geno[i]
      }
    }
}

function updateGenes(obj, callback) {
  var data_url = 'https://mygene.info/v2/query?q=' +
    'chr' + obj.chrom + ':' + obj.start + '-' + obj.end +
    '&fields=genomic_pos,symbol,name,type_of_gene'
    
  d3.json(data_url, function(err, json) {
  //d3.json('genes.json', function(err, json) {
    if (json.hits.length) {

      // Fill the data.genes array with new gene annotations.
      data.genes = new Array(json.hits.length)

      for (var i = 0; i < json.hits.length; i++) {
        if (json.hits[i].genomic_pos.length) {
          for (var j = 0; j < json.hits[i].genomic_pos.length; j++) {
            if (obj.chrom == json.hits[i].genomic_pos[j].chr) {
              data.genes[i] = {
                'id': json.hits[i]._id,
                'chrom': json.hits[i].genomic_pos[j].chr,
                'start': json.hits[i].genomic_pos[j].start,
                'end': json.hits[i].genomic_pos[j].end,
                'strand': json.hits[i].genomic_pos[j].strand,
                'symbol': json.hits[i].symbol,
                'name': json.hits[i].name,
                'type': json.hits[i].type_of_gene
              }
            }
          }
        } else {
          data.genes[i] = {
            'id': json.hits[i]._id,
            'chrom': json.hits[i].genomic_pos.chr,
            'start': json.hits[i].genomic_pos.start,
            'end': json.hits[i].genomic_pos.end,
            'strand': json.hits[i].genomic_pos.strand,
            'symbol': json.hits[i].symbol,
            'name': json.hits[i].name,
            'type': json.hits[i].type_of_gene
          }
        }
      }

      callback()
    } else {
      d3.select("#loading-text")
        .text("[mygene.info] No genes in this region: " + query)
    }
  })
}

function plotLD(target_id) {

  var margin = {top: 20, right: 20, bottom: 150, left: 40},
      width = 960 - margin.left - margin.right,
      height = 500 - margin.top - margin.bottom

  /*
  * value accessor - returns the value to encode for a given data object.
  * scale - maps value to a visual display encoding, such as a pixel position.
  * map function - maps from data value to display value
  * axis - sets up axis
  */ 

  // setup x 
  var xValue = function(d) { return d.pos }, // data -> value
      xScale = d3.scale.linear().range([0, width]), // value -> display
      xMap = function(d) { return xScale(xValue(d));}, // data -> display
      xAxis = d3.svg.axis().scale(xScale).orient("bottom"),
      xMargin = 0.1 * (d3.max(data.markers, xValue) - d3.min(data.markers, xValue));

  // setup y
  var yValue = function(d) { return d.r2 }, // data -> value
      yScale = d3.scale.linear().range([height, 0]), // value -> display
      yMap = function(d) { return yScale(yValue(d));}, // data -> display
      yAxis = d3.svg.axis().scale(yScale).orient("left");

  // setup fill color
  var cValue = function(d) { return d.chosen },
      // color = d3.scale.category20()
      color = d3.scale.ordinal()
        .range([
          '#ffffff', '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ])


  // Remove the old panel and make a new one.
  d3.select("#plotLD-svg").remove();
  var svg = d3.select(target_id).append("svg")
      .attr("id", "plotLD-svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // add the tooltip area to the webpage
  var tooltip = d3.select("body").append("div")
      .attr("class", "tooltip")
      .style("opacity", 0);

  // don't want dots overlapping axis, so add in buffer to data domain
  xScale.domain([d3.min(data.markers, xValue) - xMargin, d3.max(data.markers, xValue) + xMargin]);
  yScale.domain([d3.min(data.markers, yValue) - 0.06, d3.max(data.markers, yValue) + 0.06]);

  // x-axis
  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis)
    .append("text")
      .attr("class", "label")
      .attr("x", width)
      .attr("y", -6)
      .style("text-anchor", "end")
      .text("Position");

  // y-axis
  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis)
    .append("text")
      .attr("class", "label")
      // .attr("transform", "rotate(-90)")
      .attr("y", 6)
      .attr("x", 6)
      .attr("dy", ".71em")
      // .style("text-anchor", "end")
      .text("R²");

  var draw_legend = function() {
    var legend = svg.selectAll(".legend")
        .data(color.domain())
        .enter()
        .append("g")
        .attr("class", "legend")
        .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

    // draw legend colored rectangles
    legend.append("rect")
        .attr("x", width - 18)
        .attr("width", 18)
        .attr("height", 18)
        .style("fill", color)
        .on("click", function(d) {
            redraw(d)
        })

    // draw legend text
    legend.append("text")
        .attr("x", width - 24)
        .attr("y", 9)
        .attr("dy", ".35em")
        .style("text-anchor", "end")
        .text(function(d) { return d })
  }

  var draw_table = function() {

    // Remove the old count.
    d3.select("#marker-count").remove()
    d3.select("#main-table")
        .append("div")
        .attr("id", "marker-count")
        .text(d3.format(",")(data.markers.length) + " variants in this region")

    // Remove the old table.
    d3.select("#marker-table").remove()

    var table = d3.select("#main-table")
        .append("table")
        .attr("id", "marker-table")
        .attr("width", "100%")

    var thead = table.append("thead"),
        tbody = table.append("tbody")

    var columns = ['chrom', 'pos', 'rsid', 'ref', 'alt', 'maf', 'dist', 'r2', 'dp']

    // append the header row
    thead.append("tr")
        .selectAll("th")
        .data(columns)
        .enter()
        .append("th")
            .text(function(column) { return column })

    var sort_markers = function(a, b) {
        var x = b.r2 - a.r2
        return x == 0 ? a.dist - b.dist : x
    }

    // create a row for each object in the data
    var rows = tbody.selectAll("tr")
        .data(
            data.markers.concat().sort(sort_markers).slice(0, 15)
        )
        .enter()
        .append("tr");

    var floats = {r2:0, dp:0, maf:0}

    // create a cell in each row for each column
    var cells = rows.selectAll("td")
        .data(function(row) {
            return columns.map(function(column) {
                return {
                    column: column,
                    value: column in floats
                        ? d3.format('.2f')(row[column])
                        : row[column]
                }
            })
        })
        .enter()
        .append("td")
        .attr("font-family", "Ubuntu Mono, Monospace")
        .html(function(d) { return d.value })
    
    return table
  }

  var redraw = function(rsid) {
        tooltip.transition()
          .duration(1500)
          .style("opacity", 0);

        // Choose a reference marker that matches the rsid.
        var idx = data.markers.map(function(x) {return x.rsid}).indexOf(rsid)
        if (idx == -1) {
          idx = Math.floor(data.markers.length / 2)
          rsid = data.markers.map(function(x) {return x.rsid})[idx]
        }

        // Get the genotypes from the global 'data'.
        var geno = data.markers.map(function(x) {return x.geno})

        var ld = geno.map(function(y) {
          return compute_ld(geno[idx], y)
        })

        // // Compute D' (D prime).
        // var dp = geno.map(function(y) { 
        //   return dprime(geno[idx], y)
        // })

        // // Compute the squared Pearson correlations.
        // var r2 = geno.map(function(y) { 
        //   var r = pearson(geno[idx], y)
        //   return r * r
        // })

        // Update the global 'data'.
        for (var i = 0; i < data.markers.length; i++) {
          data.markers[i].r2 = ld[i].r2
          data.markers[i].dp = ld[i].dp
          data.markers[i].chosen = i == idx ? rsid : data.markers[i].chosen
        }

        // Transition dots to new positions.
        svg.selectAll(".dot")
            .data(data.markers)
            .transition()
            .duration(1250)
            .attr("r", function(d) { return d.chosen ? 6.0 : 3.5 })
            .attr("cx", xMap)
            .attr("cy", yMap)
            .style("fill", function(d) { return color(cValue(d)) }) 

        draw_legend()

        draw_table()
  }

  // draw dots
  svg.selectAll(".dot")
      .data(data.markers)
      .enter().append("circle")
      .attr("class", "dot")
      .attr("r", function(d) { return d.chosen ? 6.0 : 3.5 })
      .attr("cx", xMap)
      .attr("cy", yMap)
      .attr("fill", function(d) { return color(cValue(d));}) 
      .on("mouseover", function(d) {
        tooltip.transition()
          .duration(200)
          .style("opacity", .9);
        tooltip.html(
          d.rsid + "<br/>"
          + "Position: " + d3.format(",")(d.pos) + "<br/>"
          + "R²: " + d3.format(".2f")(d.r2) + "<br/>"
          + "MAF: " + d3.format(".2f")(d.maf)  
        )
        .style("left", (d3.event.pageX + 15) + "px")
        .style("top", (d3.event.pageY - 28) + "px");
      })
      .on("mouseout", function(d) {
        tooltip.transition()
          .duration(1500)
          .style("opacity", 0);
      })
      .on("click", function(d) {
        redraw(d.rsid)
      })

  draw_legend()

  draw_table()
}

function plotGenes(target_id) {

  var margin = {top: 20, right: 20, bottom: 90, left: 40},
      width = 960 - margin.left - margin.right,
      height = 500 - margin.top - margin.bottom

  /*
  * value accessor - returns the value to encode for a given data object.
  * scale - maps value to a visual display encoding, such as a pixel position.
  * map function - maps from data value to display value
  * axis - sets up axis
  */ 

  // setup x 
  var xValue = function(d) { return d.pos }, // data -> value
      xScale = d3.scale.linear().range([0, width]), // value -> display
      xMap = function(d) { return xScale(xValue(d));}, // data -> display
      xAxis = d3.svg.axis().scale(xScale).orient("bottom"),
      xMargin = 0.1 * (d3.max(data.markers, xValue) - d3.min(data.markers, xValue));

  // setup y
  var yValue = function(d) { return d.r2 }, // data -> value
      yScale = d3.scale.linear().range([height, 0]), // value -> display
      yMap = function(d) { return yScale(yValue(d));}, // data -> display
      yAxis = d3.svg.axis().scale(yScale).orient("left");

  // setup fill color
  var cValue = function(d) { return d.type },
      color = d3.scale.category20()
      // color = d3.scale.ordinal()
      //   .range([
      //     '#ff7f0e', '#1f77b4', '#2ca02c', '#d62728', '#9467bd',
      //     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
      //   ])


  // Remove the old panel and make a new one.
  var svg = d3.select("#plotLD-svg")

  // don't want dots overlapping axis, so add in buffer to data domain
  xScale.domain([d3.min(data.markers, xValue) - xMargin, d3.max(data.markers, xValue) + xMargin]);
  yScale.domain([d3.min(data.markers, yValue) - 0.06, d3.max(data.markers, yValue) + 0.06]);

  // TODO: Use getComputedTextLength instead of "20 * length"
  // var ty = svg.append("text")
  //     .text(msg);
  //     // could also use getBBox for width and height
  //     var wid = ty.node().getComputedTextLength();
  
  // Use an interval tree to place them in non-overlapping positions.
  var it = [new IntervalTree()],
      n_levels = 5
  for (var i = 0; i < data.genes.length; i++) {
    // var s = data.genes[i].start - 5e3,
    //     e = data.genes[i].end + 5e3
    var s = xScale(data.genes[i].start)
    var e = Math.max(xScale(data.genes[i].end), s + 20 * data.genes[i].symbol.length)
    for (var level = 0; level < n_levels; level++) {
      if (!it[level]) {
        it.push(new IntervalTree())
      }
      if (!it[level].intersects([s,e])) {
        it[level].add([s,e])
        data.genes[i].level = level
        break
      }
    }
  }

  svg.selectAll(".gene")
    .data(data.genes)
    .enter()
    .append("rect")
    .attr("class", "gene")
    .attr("x", function(d) { return xScale(d.start) })
    .attr("width", function(d) { return Math.abs(xScale(d.end) - xScale(d.start)) })
    //.attr("y", yScale(-0.075))
    .attr("y", function(d) { return height + 30 * d.level - 2 })
    .attr("height", 6)
    .attr("fill", function(d) { return color(cValue(d)) })
    .on("click", function(d) {
      console.log(d)
    })

  svg.selectAll(".gene-label")
    .data(data.genes)
    .enter()
    .append("text")
    .attr("class", "gene-label")
    .attr("x", function(d) { return Math.max(0, xScale(d.start)) })
    .attr("y", function(d) { return height + 30 * d.level + 20 })
    .text(function(d) { return d.symbol })

  //<rect x="0" y="0" width="50" height="50" fill="green" />
}

// Pure functions ------------------------------------------------------------

// Return a URL to the 1000 Genomes data hosted at Broad Institute.
function dataURL(chrom, start, end, limit) {
  //return 'https://data.broadinstitute.org/srlab/BEAGLE/1000_Genomes_phase3_v5a/chr' +
  //  chrom + '.1kg.phase3.v5a.vcf.gz?' + encodeData({
  //    format: 'text',
  //    limit: limit,
  //    region: `${chrom}:${start}-${end}`
  //  })
  // BEAGLE webserver
  var host_url = 'http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/'
  var filename = `chr${chrom}.1kg.phase3.v5a.vcf.gz`
  // Amazon AWS
  // var host_url = 'http://s3.amazonaws.com/1000genomes/release/20101123/interim_phase1_release/'
  // var filename = `ALL.chr${chrom}.phase1.projectConsensus.genotypes.vcf.gz`
  var region = `${chrom}:${start}-${end}`
  var tabix_url = `https://tabix.iobio.io/?cmd=-h%20%27${host_url}${filename}%27%20${region}`
  return tabix_url
}

// Return the genotypes in a VCF row:
//      genotypes([ 
//          ['6', '528156', 'rs141711206', 'A', 'T', '.', 'PASS', '.', 'GT', '0|1 1|1 0|0...']
//      ])
// as a numeric vector:
//      [0, 1, 1, 1, 0, 0, ...]
function genotypes(vcf_row) {
    return vcf_row.slice(10,-1).join("|").split("|").map(function(x) {
        return +x
    })
}

// Subtract the mean from a vector of numbers.
function center(xs) {
    var avg = d3.sum(xs) / xs.length
    return xs.map(function(x) { return x - avg })
}

// Return the dot product of two vectors.
function dot(x, y) {
  var sum = 0
  for (var i = 0; i < x.length; i++) {
      sum += x[i] * y[i]
  }
  return sum
}

// x is a vector of genotypes
// y is a vector of genotypes
// Encoding is 0 for reference, 1 for alternate.
// Each person has 2 alleles.
// For example, here is a vector for 3 people:
//    [0, 0, 0, 1, 1, 1]
// The first person is homozygous for the reference allele.
function compute_ld(x, y) {
  var retval = {r2: 0, dp: 0}
  // These computations are used for both r2 and dp.
  var px = d3.sum(x) / x.length,
      py = d3.sum(y) / y.length,
      pxy = 0
  for (var i = 0; i < x.length; i++) {
    pxy += x[i] && y[i]
  }
  pxy /= x.length
  var d = pxy - px * py,
      dmax = 0
  if (d < 0) {
    dmax = Math.min(px * py, (1 - px) * (1 - py))
  } else {
    dmax = Math.min(px * (1 - py), (1 - px) * py)
  }
  // D prime.
  if (dmax != 0) {
    retval.dp = Math.abs(d / dmax)
  } else {
    retval.dp = 0
  }
  // Squared Pearson correlation coefficient.
  retval.r2 = Math.pow(d / Math.sqrt(px * (1 - px) * py * (1 - py)), 2)
  return retval
}

// // Return the Pearson correlation coefficient of two vectors.
// // x = [1.2,1.3,1,0,0.1,0]
// // y = [1.5,0.4,1,0,1.1,0]
// // function pearson(x, y) {
// //   var num = dot(x, y),
// //       den = Math.sqrt(dot(x, x) * dot(y, y))
// //   return den ? num / den : 1
// // }
// function pearson(x, y) {
//   var px = d3.sum(x) / x.length,
//       py = d3.sum(y) / y.length,
//       pxy = 0
//   for (var i = 0; i < x.length; i++) {
//     pxy += x[i] && y[i]
//   }
//   pxy /= x.length
//   var d = pxy - px * py
//   return d / Math.sqrt(px * (1 - px) * py * (1 - py))
// }
// 
// // Return the D prime of two bi-allelic markers.
// // x = [1,1,1,0,0,0]
// // y = [1,0,1,0,1,0]
// function dprime(x, y) {
//   // var x = [1,1,1,0,0,0],
//   //     y = [1,0,1,0,1,0]
//   var px = d3.sum(x) / x.length,
//       py = d3.sum(y) / y.length,
//       pxy = 0
//   for (var i = 0; i < x.length; i++) {
//     pxy += x[i] && y[i]
//   }
//   pxy /= x.length
//   var d = pxy - px * py,
//       dmax = 0
//   if (d < 0) {
//     dmax = Math.min(px * py, (1 - px) * (1 - py))
//   } else {
//     dmax = Math.min(px * (1 - py), (1 - px) * py)
//   }
//   if (dmax != 0) {
//     return d / dmax
//   }
//   return 0
// }

// Encode an object with keys and values into a string for an URL.
// Calls encodeURIComponent() on each key and value.
// E.g.
//      encodeDate({limit:-1,rsid:'rs123'})
//      "limit=-1&rsid=123"
function encodeData(obj) {
  return Object.keys(obj).map(function(key) {
    return [key, obj[key]].map(encodeURIComponent).join("=");
  }).join("&");
}   

// Parse a region "3:1,230-4,565" into {chrom: '3', start: 1230, end: 4565}
function parseRegion(region) {
  if (region.indexOf(':') > -1) {
    var chrom = region.split(':')[0]
    var pos = region.split(':')[1].replace(/,/g, '')
    var start = pos
    var end = pos
    if (pos.indexOf('-') > -1) {
      start = +pos.split('-')[0]
      end = +pos.split('-')[1]
    }
    if (start == end) {
      start = +start - 50000
      end = +end + 50000
    }
    return {
      chrom: chrom,
      start: start,
      end: end
    }
  }
  return region
}

// Create a region "3:1,230-4,565" with toRegion("3", 1230, 4565)
function toRegion(chrom, start, end) {
  var cf = d3.format(",")
  return chrom + ':' + cf(start) + '-' + cf(end)
}

function shuffle(array) {
  var counter = array.length, temp, index;

  // While there are elements in the array
  while (counter > 0) {
    // Pick a random index
    index = Math.floor(Math.random() * counter);

    // Decrease counter by 1
    counter--;

    // And swap the last element with it
    temp = array[counter];
    array[counter] = array[index];
    array[index] = temp;
  }

  return array;
}

/*
function roundNumber(x, digits) {
  digits = typeof digis !== 'undefined' ? digits : 2;
  if (x === 0) {
    return 0;
  }
  if (x < 0.01) {
    return x.toExponential(digits);
  }
  return x.toPrecision(digits + 1);
};

function seq(n) {
  retval = []
  for (var i = 0; i < n; i++) {
    retval.push(i);
  }
  return retval;
}
*/
