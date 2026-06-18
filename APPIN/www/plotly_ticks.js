/* =============================================================================
   2D NMR Analyst - Plotly Tick Management
   Author: Julien Guibert
   
   This script handles synchronous tick updates during zoom operations
   on the interactive NMR spectrum plot.
   ============================================================================= */

/**
 * Generate 'clean' tick values for a given range
 * Uses NMR convention (inverted sign for display)
 * 
 * @param {number} min - Minimum value of the range
 * @param {number} max - Maximum value of the range
 * @param {number} targetCount - Target number of ticks
 * @returns {Object} Object with tickvals and ticktext arrays
 */
function generateNiceTicks(min, max, targetCount) {
  var range = max - min;
  if (range <= 0) return { tickvals: [], ticktext: [] };
  
  // Nice step values
  var niceSteps = [0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50, 100];
  var roughStep = range / targetCount;
  
  // Find the closest nice step
  var step = niceSteps[0];
  for (var i = 0; i < niceSteps.length; i++) {
    if (niceSteps[i] >= roughStep) {
      step = niceSteps[i];
      break;
    }
    step = niceSteps[i];
  }
  
  // Calculate ticks
  var startTick = Math.ceil(min / step) * step;
  var tickvals = [];
  var ticktext = [];
  
  // Determine number of decimals
  var decimals = 0;
  if (step < 1) decimals = step < 0.1 ? 2 : 1;
  
  for (var t = startTick; t <= max; t += step) {
    tickvals.push(t);
    // Invert sign for display (NMR convention)
    ticktext.push((-t).toFixed(decimals));
  }
  
  return { tickvals: tickvals, ticktext: ticktext };
}

/**
 * Main function to update ticks on zoom
 * 
 * @param {HTMLElement} gd - The Plotly graph div element
 */
function updateTicksOnZoom(gd) {
  if (!gd || !gd.layout) return;
  
  var xaxis = gd.layout.xaxis || {};
  var yaxis = gd.layout.yaxis || {};
  
  var xRange = xaxis.range;
  var yRange = yaxis.range;
  
  if (!xRange || !yRange) return;
  
  var xTicks = generateNiceTicks(Math.min(xRange[0], xRange[1]), Math.max(xRange[0], xRange[1]), 10);
  var yTicks = generateNiceTicks(Math.min(yRange[0], yRange[1]), Math.max(yRange[0], yRange[1]), 10);
  
  // Synchronous update without triggering relayout event
  Plotly.relayout(gd, {
    'xaxis.tickmode': 'array',
    'xaxis.tickvals': xTicks.tickvals,
    'xaxis.ticktext': xTicks.ticktext,
    'xaxis.showticklabels': true,
    'xaxis.ticks': 'outside',
    'yaxis.tickmode': 'array',
    'yaxis.tickvals': yTicks.tickvals,
    'yaxis.ticktext': yTicks.ticktext,
    'yaxis.showticklabels': true,
    'yaxis.ticks': 'outside'
  });
}

/**
 * Hide ticks temporarily (used during unzoom/reset)
 * 
 * @param {HTMLElement} gd - The Plotly graph div element
 */
function hideTicksTemporarily(gd) {
  Plotly.relayout(gd, {
    'xaxis.showticklabels': false,
    'xaxis.ticks': '',
    'yaxis.showticklabels': false,
    'yaxis.ticks': ''
  });
}

/**
 * Observer to detect when the plot is created/updated
 * Attaches zoom event listeners to the interactive plot
 */
$(document).on('shiny:value', function(event) {
  if (event.name === 'interactivePlot') {
    setTimeout(function() {
      var gd = document.getElementById('interactivePlot');
      if (gd && gd.on) {
        // Listen for zoom events
        gd.on('plotly_relayout', function(eventData) {
          if (!eventData) return;
          
          // Ignore our own tick updates
          if (eventData['xaxis.tickvals'] !== undefined) return;
          if (eventData['xaxis.showticklabels'] !== undefined && 
              eventData['xaxis.range[0]'] === undefined &&
              eventData['xaxis.autorange'] === undefined) return;
          
          // Detect if it's an autoscale (double-click, reset)
          var isAutoscale = eventData['xaxis.autorange'] !== undefined || 
                            eventData['yaxis.autorange'] !== undefined ||
                            eventData['autosize'] !== undefined;
          
          if (isAutoscale) {
            // Hide ticks immediately during autoscale
            hideTicksTemporarily(gd);
            // Wait for plotly to calculate new ranges then update
            setTimeout(function() { updateTicksOnZoom(gd); }, 80);
          } else if (eventData['xaxis.range[0]'] !== undefined || 
                     eventData['yaxis.range[0]'] !== undefined) {
            // Manual zoom - unchanged behavior
            setTimeout(function() { updateTicksOnZoom(gd); }, 10);
          }
        });
      }
    }, 500);
  }
});
