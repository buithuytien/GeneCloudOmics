[ {"selector": "node", "css": {
      "shape": "ellipse",
      "text-valign":"center",
      "text-halign":"center",
      "content": "data(id)",
      "background-color": "#FFFFFF",
      "border-color": "black","border-width":"1px",
      "width":  "40px",
      "height": "40px",
      "font-size":"14px"}},

   {"selector":"node:selected", "css": {
       "text-valign":"center",
       "text-halign":"center",
       "border-color": "black",
       "overlay-opacity": 0.2,
       "overlay-color": "gray"
       }},

    {"selector": "node[[degree >= 5]]", "css": {
        "width": "60px",
        "height": "60px"
    }},

    {"selector": "node[[degree >= 8]]", "css": {
        "width": "70px",
        "height": "70px"
    }},

    {"selector": "node[[degree >= 12]]", "css": {
        "width": "80px",
        "height": "80px"
    }},

   {"selector": "edge", "css": {
       "opacity": 0.5,
       "curve-style": "bezier"
       }},

   {"selector":"edge:selected", "css": {
       "overlay-opacity": 0.2,
       "overlay-color": "red"
       }}
]
