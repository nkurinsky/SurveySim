FUNCTION COLOR24, number

   ; This FUNCTION accepts a [red, green, blue] triple that
   ; describes a particular color and returns a 24-bit long
   ; integer that is equivalent to that color. The color is
   ; described in terms of a hexidecimal number (e.g., FF206A)
   ; where the left two digits represent the blue color, the
   ; middle two digits represent the green color, and the right
   ; two digits represent the red color.
   ;
   ; The triple can be either a row or column vector of 3 elements.

ON_ERROR, 1

IF N_ELEMENTS(number) NE 3 THEN $
   MESSAGE, 'Augument must be a three-element vector.'

IF MAX(number) GT 255 OR MIN(number) LT 0 THEN $
   MESSAGE, 'Argument values must be in range of 0-255'

base16 = [[1L, 16L], [256L, 4096L], [65536L, 1048576L]]

num24bit = 0L

FOR j=0,2 DO num24bit = num24bit + ((number(j) MOD 16) * base16(0,j)) + $
   (Fix(number(j)/16) * base16(1,j))

RETURN, num24bit
END ; ************************  of COLOR24  ******************************


FUNCTION GETCOLOR, thisColor, index, TRUE=truecolor, $
   NAMES=colornames, LOAD=load, START=start, INDEXED=indexedcolor

   ; Set up the color vectors.

   names =  ['White']
   rvalue = [ 255]
   gvalue = [ 255]
   bvalue = [ 255]
   names  = [ names,        'Snow',     'Ivory','Light Yellow',   'Cornsilk',      'Beige',   'Seashell' ]
   rvalue = [ rvalue,          255,          255,          255,          255,          245,          255 ]
   gvalue = [ gvalue,          250,          255,          255,          248,          245,          245 ]
   bvalue = [ bvalue,          250,          240,          224,          220,          220,          238 ]
   names  = [ names,       'Linen','Antique White',    'Papaya',     'Almond',     'Bisque',  'Moccasin' ]
   rvalue = [ rvalue,          250,          250,          255,          255,          255,          255 ]
   gvalue = [ gvalue,          240,          235,          239,          235,          228,          228 ]
   bvalue = [ bvalue,          230,          215,          213,          205,          196,          181 ]
   names  = [ names,       'Wheat',  'Burlywood',        'Tan', 'Light Gray',   'Lavender','Medium Gray' ]
   rvalue = [ rvalue,          245,          222,          210,          230,          230,          210 ]
   gvalue = [ gvalue,          222,          184,          180,          230,          230,          210 ]
   bvalue = [ bvalue,          179,          135,          140,          230,          250,          210 ]
   names  = [ names,        'Gray', 'Slate Gray',  'Dark Gray',   'Charcoal',      'Black', 'Light Cyan' ]
   rvalue = [ rvalue,          190,          112,          110,           70,            0,          224 ]
   gvalue = [ gvalue,          190,          128,          110,           70,            0,          255 ]
   bvalue = [ bvalue,          190,          144,          110,           70,            0,          255 ]
   names  = [ names, 'Powder Blue',   'Sky Blue', 'Steel Blue','Dodger Blue', 'Royal Blue',       'Blue' ]
   rvalue = [ rvalue,          176,          135,           70,           30,           65,            0 ]
   gvalue = [ gvalue,          224,          206,          130,          144,          105,            0 ]
   bvalue = [ bvalue,          230,          235,          180,          255,          225,          255 ]
   names  = [ names,        'Navy',   'Honeydew', 'Pale Green','Aquamarine','Spring Green',       'Cyan' ]
   rvalue = [ rvalue,            0,          240,          152,          127,            0,            0 ]
   gvalue = [ gvalue,            0,          255,          251,          255,          250,          255 ]
   bvalue = [ bvalue,          128,          240,          152,          212,          154,          255 ]
   names  = [ names,   'Turquoise', 'Sea Green','Forest Green','Green Yellow','Chartreuse', 'Lawn Green' ]
   rvalue = [ rvalue,           64,           46,           34,          173,          127,          124 ]
   gvalue = [ gvalue,          224,          139,          139,          255,          255,          252 ]
   bvalue = [ bvalue,          208,           87,           34,           47,            0,            0 ]
   names  = [ names,       'Green', 'Lime Green', 'Olive Drab',     'Olive','Dark Green','Pale Goldenrod']
   rvalue = [ rvalue,            0,           50,          107,           85,            0,          238 ]
   gvalue = [ gvalue,          255,          205,          142,          107,          100,          232 ]
   bvalue = [ bvalue,            0,           50,           35,           47,            0,          170 ]
   names  = [ names,       'Khaki', 'Dark Khaki',     'Yellow',       'Gold','Goldenrod','Dark Goldenrod']
   rvalue = [ rvalue,          240,          189,          255,          255,          218,          184 ]
   gvalue = [ gvalue,          230,          183,          255,          215,          165,          134 ]
   bvalue = [ bvalue,          140,          107,            0,            0,           32,           11 ]
   names  = [ names,'Saddle Brown',       'Rose',       'Pink', 'Rosy Brown','Sandy Brown',       'Peru' ]
   rvalue = [ rvalue,          139,          255,          255,          188,          244,          205 ]
   gvalue = [ gvalue,           69,          228,          192,          143,          164,          133 ]
   bvalue = [ bvalue,           19,          225,          203,          143,           96,           63 ]
   names  = [ names,  'Indian Red',  'Chocolate',     'Sienna','Dark Salmon',    'Salmon','Light Salmon' ]
   rvalue = [ rvalue,          205,          210,          160,          233,          250,          255 ]
   gvalue = [ gvalue,           92,          105,           82,          150,          128,          160 ]
   bvalue = [ bvalue,           92,           30,           45,          122,          114,          122 ]
   names  = [ names,      'Orange',      'Coral', 'Light Coral',  'Firebrick',      'Brown',  'Hot Pink' ]
   rvalue = [ rvalue,          255,          255,          240,          178,          165,          255 ]
   gvalue = [ gvalue,          165,          127,          128,           34,           42,          105 ]
   bvalue = [ bvalue,            0,           80,          128,           34,           42,          180 ]
   names  = [ names,   'Deep Pink',    'Magenta',     'Tomato', 'Orange Red',        'Red', 'Violet Red' ]
   rvalue = [ rvalue,          255,          255,          255,          255,          255,          208 ]
   gvalue = [ gvalue,           20,            0,           99,           69,            0,           32 ]
   bvalue = [ bvalue,          147,          255,           71,            0,            0,          144 ]
   names  = [ names,      'Maroon',    'Thistle',       'Plum',     'Violet',    'Orchid','Medium Orchid']
   rvalue = [ rvalue,          176,          216,          221,          238,          218,          186 ]
   gvalue = [ gvalue,           48,          191,          160,          130,          112,           85 ]
   bvalue = [ bvalue,           96,          216,          221,          238,          214,          211 ]
   names  = [ names, 'Dark Orchid','Blue Violet',     'Purple' ]
   rvalue = [ rvalue,          153,          138,          160 ]
   gvalue = [ gvalue,           50,           43,           32 ]
   bvalue = [ bvalue,          204,          226,          240 ]

   ; Did the user ask for a specific color? If not, return
   ; all the colors. If the user asked for a specific color,
   ; find out if a 24-bit value is required. Return to caller
   ; if an error occurs.

ON_Error, 2
ncolors = N_Elements(names)

np = N_Params()
IF N_Elements(start) EQ 0 THEN start = !D.TABLE_SIZE - ncolors - 1 ELSE start = start < (!D.TABLE_SIZE - ncolors - 1)

   ; User ask for the color names?

IF Keyword_Set(colornames) THEN RETURN, Reform(names, 1, N_Elements(names)) $
    ELSE names = StrUpCase(StrCompress(StrTrim(names,2), /Remove_All))

   ; If no positional parameter, return all colors.

IF np EQ 0 THEN BEGIN

   ; Did the user want a 24-bit value? If so, call COLOR24.

   IF Keyword_Set(trueColor) THEN BEGIN
      returnColor = LonArr(ncolors)
      FOR j=0,ncolors-1 DO returnColor[j] = Color24([rvalue[j], gvalue[j], bvalue[j]])

         ; If LOAD keyword set, return a color structure.

      IF Keyword_Set(load) THEN BEGIN
         returnValue = Create_Struct('white', returnColor[0])
         FOR j=1,ncolors-1 DO returnValue = Create_Struct(returnValue, names[j], returnColor[j])
         returnColor = returnValue
      ENDIF

      RETURN, returnColor
   ENDIF

   ; If color decomposition is ON and INDEXED is not set, return 24-bit values.

   IF Float(!Version.Release) GE 5.2 THEN BEGIN
      IF (!D.Name EQ 'X' OR !D.Name EQ 'WIN' OR !D.Name EQ 'MAC') THEN BEGIN
         Device, Get_Decomposed=decomposedState
      ENDIF ELSE decomposedState = 0
      IF Keyword_Set(indexedcolor) THEN decomposedState = 0
      IF decomposedState EQ 1 THEN BEGIN
         returnColor = LonArr(ncolors)
         FOR j=0,ncolors-1 DO returnColor[j] = Color24([rvalue[j], gvalue[j], bvalue[j]])
         IF Keyword_Set(load) THEN BEGIN
            returnValue = Create_Struct('white', returnColor[0])
            FOR j=1,ncolors-1 DO returnValue = Create_Struct(returnValue, names[j], returnColor[j])
            RETURN, returnValue
         ENDIF
         RETURN, returnColor
      ENDIF

      IF Keyword_Set(load) THEN BEGIN
         TVLCT, Reform([rvalue, gvalue, bvalue], ncolors, 3), start
         returnValue = Create_Struct('white', start)
         FOR j=1,ncolors-1 DO returnValue = Create_Struct(returnValue, names[j], start+j)
         RETURN, returnValue
      ENDIF

      returnColor = REFORM([rvalue, gvalue, bvalue], ncolors, 3)
      RETURN, returnColor

   ENDIF

   IF Keyword_Set(load) THEN BEGIN
      TVLCT, Reform([rvalue, gvalue, bvalue], ncolors, 3), start
      returnValue = Create_Struct('white', start)
      FOR j=1,ncolors-1 DO returnValue = Create_Struct(returnValue, names[j], start+j)
      RETURN, returnValue
   ENDIF

   returnColor = REFORM([rvalue, gvalue, bvalue], ncolors, 3)
   RETURN, returnColor

ENDIF

   ; Make sure the color parameter is an uppercase string.

varInfo = SIZE(thisColor)
IF varInfo(varInfo(0) + 1) NE 7 THEN $
   MESSAGE, 'The color name must be a string.'
thisColor = StrUpCase(StrCompress(StrTrim(thisColor,2), /Remove_All))

   ; Check synonyms of colors.

IF StrUpCase(thisColor) EQ 'GREY' THEN thisColor = 'GRAY'
IF StrUpCase(thisColor) EQ 'AQUA' THEN thisColor = 'AQUAMARINE'
IF StrUpCase(thisColor) EQ 'SKYBLUE' THEN thisColor = 'SKY BLUE'
IF StrUpCase(thisColor) EQ 'LIGHTGREY' THEN thisColor = 'LIGHTGRAY'
IF StrUpCase(thisColor) EQ 'MEDIUMGREY' THEN thisColor = 'MEDIUMGRAY'
IF StrUpCase(thisColor) EQ 'SLATEGREY' THEN thisColor = 'SLATEGRAY'
IF StrUpCase(thisColor) EQ 'DARKGREY' THEN thisColor = 'DARKGRAY'
IF StrUpCase(thisColor) EQ 'SKY' THEN thisColor = 'SKY BLUE'
IF StrUpCase(thisColor) EQ 'NAVY BLUE' THEN thisColor = 'NAVY'
IF StrUpCase(thisColor) EQ 'NAVYBLUE' THEN thisColor = 'NAVY'


   ; Get the color triple for this color.

colorIndex = WHERE(names EQ thisColor)

   ; If you can't find it. Issue an informational message,
   ; set the index to a YELLOW color, and continue.

IF colorIndex(0) LT 0 THEN BEGIN
   MESSAGE, "Can't find color " + thisColor + ". Returning " + StrUpCase(names[0]) + ".", /INFORMATIONAL
   thisColor = names[0]
   colorIndex = 0
ENDIF

   ; Get the color triple.

r = rvalue(colorIndex)
g = gvalue(colorIndex)
b = bvalue(colorIndex)
returnColor = REFORM([r, g, b], 1, 3)

   ; Did the user want a 24-bit value? If so, call COLOR24.

IF KEYWORD_SET(trueColor) THEN BEGIN
   returnColor = COLOR24(returnColor)
   RETURN, returnColor[0]
ENDIF

   ; If color decomposition is ON and INDEXED is OFF,, return 24-bit value.

IF Float(!Version.Release) GE 5.2 THEN BEGIN

   IF (!D.Name EQ 'X' OR !D.Name EQ 'WIN' OR !D.Name EQ 'MAC') THEN BEGIN
      Device, Get_Decomposed=decomposedState
   ENDIF ELSE decomposedState = 0
   IF Keyword_Set(indexedcolor) THEN decomposedState = 0

   IF decomposedState EQ 1 THEN BEGIN

         ; Before you change return color, load index if requested.

      IF N_Elements(index) NE 0 THEN BEGIN
         index = 0 > index < (!D.Table_Size-1)
         TVLCT, returnColor, index
      ENDIF

      returnColor = COLOR24(returnColor)
      RETURN, returnColor[0]
   ENDIF
ENDIF

   ; Did the user specify a color index? If so, load it.

IF N_Elements(index) NE 0 THEN BEGIN
   index = 0 > index < (!D.Table_Size-1)
   TVLCT, returnColor, index
   returnColor = index
   RETURN, returnColor[0]
ENDIF

   ; Did the user specify INDEXED color? If so, load it.

IF Keyword_Set(indexedColor) THEN BEGIN
   TVLCT, returnColor, !P.Color
   returnColor = !P.Color < (!D.Table_Size -1)
   RETURN, returnColor[0]
ENDIF

RETURN, returnColor
END
