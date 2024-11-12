#' Einlesealgorithmus für klinisches XML-Format des Zentrum für Krebsregisterdaten (ZfKD) am Robert-Koch-Institut (RKI)
#' Herausgeber: Institut für digitale Gesundheitsdaten Rheinland-Pfalz gGmbH
#' Autor: Benjamin Rieger (`rieger@idg-rlp.de`)
#' Stand: November 2024
#' Repository: https://github.com/IDGRLP/ParseXML


# Libraries ####
library(dplyr)
library(data.table)
library(rlang)
library(purrr)
library(hablar)
library(xml2)
library(gdata)


# Functions ####
# Gib redundanten Elementen in einem Vektor eine Laufnummer
enumEntsInVec = function(.x, .sep = ':', .append = TRUE, .single.tag = '0'){
  .enum.vec = character(length = length(.x))
  
  for(i in unique(.x[!is.na(.x)])){
    .mask = .x == i & !is.na(.x)
    .range = sum(.mask)
    
    if(.range > 1)
      .enum.vec[.mask] = as.character(1:.range)
    else
      .enum.vec[.mask] = as.character(.single.tag)
  }
  
  if(.append)
    paste0(.x, .sep, .enum.vec)
  else
    .enum.vec
}



# Lade xml
xml <- xml2::read_xml('<< Enter here XML-File.xml >>')
# holen uns den Teil der Subliste die uns interessiert
lxml = xml2::as_list(xml2::xml_children(xml)[[3]])


# Rekursionsfunktion
# .actual.id = aktueller Knotenname in dem die Funktion sich befindet
# .enum.id = Laufnummer des aktuellen Knotens (1-x, Falls ID oefters vorkommt, '0' bei Einzelknoten)
# .add.id.to.path = FLAG ob Knoten in Pfad aufgenommen werden soll
# .path = absoluter Pfad, wird stetig fuer jede Rekursion weiter gebaut
# .level = Level der Rekursionstiefe, +1 fuer jeden weiteren Durchlauf (Debug)
recList = function(.list, .actual.id = '', .enum.id = '', .add.id.to.path = FALSE, .path = NULL, .level = 1){
  # alle Attributswerte des aktuellen Listenelements (named vector)
  .atts = attributes(.list)
  # alle Attributsnamen
  .atts.ids = names(.atts)
  
  # Falls es weitere (Sub-)Knoten gibt (Attribut 'names' ist gesetzt und gefuellt; beim Endknoten leer)...
  if('names' %in% .atts.ids){
    # Falls der Knoten in den Pfad eingebracht werden soll (da es weitere Subknoten gibt (-> 'names), wird der nummerierte Knotennamen (.enum.id) eingebracht)
    if(.add.id.to.path)
      .path = c(.path, paste0(.actual.id, ':', .enum.id))
    
    # Attribut 'names' wird entfernt...
    .atts.ids = .atts.ids[!.atts.ids %in% c('names')]
    
    # Falls es weitere Attribute gibt (Attribute koennen wie Knoten in den Pfad eingepflegt werden)...
    if(length(.atts.ids) > 0){
      # ...baue key:value vektor auf
      .key.value = sort(paste0(.atts.ids, ':', .atts[.atts.ids]))
      # ...und pflege in den Pfad ein
      .path = c(.path, .key.value)
    }
    
    # hole alle Sub-Knoten-IDs aus Attribut 'names'
    .ids = names(.list)
    
    # Erzeuge Nummerierungsvektor fuer Sub-Knoten-IDs
    # Diese dienen den Knoten als eindeutiger Identifizierungswert fuer die spaetere Ueberfuehrung in die Tabellen
    .ids.enum = enumEntsInVec(.ids, .append = FALSE, .single.tag = '0')
  
    # setze '.add.id.to.path' = TRUE, ab jetzt werden alle Knoten in den Pfad geschrieben (verhindert das Setzen der leeren Knoten-ID des ersten Aufrufs)
    .add.id.to.path = TRUE
    
    # iteriere ueber Positionen der Knotennamen
    for(.pos in seq_along(.ids))
      # rufe Rekursion mit Sub-Knoten, dessen ID, Nummerierung, und Flags auf
      recList(.list[[.pos]], .ids[[.pos]], .ids.enum[[.pos]], .add.id.to.path, .path, .level + 1)
    # Falls wir am Ende eines Knotens sind ('names'-Attribut war NICHT gesetzt)...
  } else{
    # letztes Element (Value) als Liste vorhanden
    # ...setze die (Value-)Liste zu einem String zusammen (Kommasepariert)
    .value = paste0(sort(unlist(.list)), collapse = ', ')

    # ...falls der Endknoten auf diesem Level unique war (Nummerierungstag == '0')
    if(.enum.id %in% '0')
      # ...Nimm die einfache ID fuer den key:value string
      .key.value = paste0(.actual.id, ':', .value)
    # ...falls es den Endknoten oefters gab...
    else
      # ...haenge die Nummerierung feste am Namen an
      .key.value = paste0(paste0(.actual.id, '_', .enum.id), ':', .value)
    
    # ...setze den Pfad-String zusammen (Pfad bis hierher)
    .path.string = paste0(paste0(.path, collapse = ' # '), ' # ', .key.value)

    # ...speichere Pfad in LookUp
    LookUp[[.path.string]] = 1L
    
    # Falls es weitere Attribute in den LeafNodes gibt (Attribut 'names' war im else-Block nicht dabei)...
    if(length(.atts.ids) > 0){
      # ...baue key:value Vektor auf, ergaenze keys mit aktueller Knoten_ID (um zugehoerigkeit zu zeigen; -> <KnotenID>_<key>:<value>)
      .key.value = paste0(.actual.id, '_', .atts.ids, ':', .atts[.atts.ids])

      # iteriere ueber .key.value...
      walk(.key.value, ~{
        # ...baue Pfad auf (key:value wurde angepasst)
        .path.string = paste0(paste0(.path, collapse = ' # '), ' # ', .x)
        # ...& speichere Pfad
        LookUp[[.path.string]] = 1L
      })
    }
  }
}

# Erstelle LookUp
LookUp = new.env(hash = TRUE)
# Laufe Rekursion
recList(lxml)
# ueberpruefe Ergebnis #
# print(sort(ls(LookUp)))



# Knoten, welche als Schluessel zur Identifikation der Werte dienen, je mehr Schluesselknoten angegeben, desto mehr Tabellen koennen entstehen #
keyNodes = c('Patient_ID', 'Tumor_ID', 'Primaerdiagnose', 'OP', 'ST', 'SYST', 'Folgeereignis', 'Weitere_Klassifikation',
             'Modul_Darm', 'Modul_Malignes_Melanom', 'Modul_Mamma', 'Modul_Prostata')

# Knoten welche aus dem Pfad entfernt werden koennen
ignoreNodes = c('Patient', 'Patienten_Stammdaten', 'Vitalstatus', 'Todesursachen', 'Tumor',
                'Menge_Tumor', 'Menge_SYST', 'Menge_ST', 'Menge_OP', 'Menge_Folgeereignis',
                'Menge_Weitere_Klassifikation', 'Menge_FM', 'Menge_Substanz', 'Menge_Bestrahlung', 'Menge_OPS', 'Menge_Weitere_Todesursachen',
                'Zielgebiet', 'Applikationsart')


# Erstelle LookUp
Aggregate.LookUp = new.env(hash = TRUE)

# iteriere ueber Pfade...
walk(sort(ls(LookUp)), ~{
  # Knoten sind durch ' # ' getrennt und die Werte haengen mit ':' an den jeweiligen Knoten an (key:value # key:value # key:value # .....)
  
  # zerlege den Keystring in key-value-Liste
  .key.value.vec = unlist(strsplit(.x, ' # ')) %>% strsplit(':')
  # erzeuge named vector mit keys als Namen & values als Werten
  .frame = map_chr(.key.value.vec, ~.x[[2]])
  names(.frame) = map_chr(.key.value.vec, ~.x[[1]])
  
  # erhalte Wert der Leaf-Nodes (der Wert welcher in jedem Pfad am Ende steht)
  .value = .frame[length(.frame)]

  # entferne unerwuenschte Knoten (-> ignoreNodes)
  .frame = .frame[!names(.frame) %in% ignoreNodes]
  # erzeuge Namensvektor
  .frame.id = names(.frame)

  # Markiere Schluesselknoten in bool-vektor
  .keynode.mask = .frame.id %in% keyNodes
  # berechne Position des am weitesten im Pfad liegenden Schluesselknotens
  .knm.max = max_(which(.keynode.mask))
  
  # Falls letzter Schluesselknoten nicht letzter Knoten ist...
  if(.knm.max < length(.frame)){
    # erzeuge aus Schluesselknoten wieder key:value string (separiert mit '#' bei mehreren key:value Paaren)
    .tab.key = paste0(.frame.id[.keynode.mask], ':', .frame[.keynode.mask], collapse = ' # ')
    
    # Erzeuge flattenNode (alle Knotennamen nach dem letzten Schluessel werden konkateniert)
    # falls es nur noch um den LeafNode geht...
    if(.knm.max + 1 == length(.frame)){
      # ...setze flattenNode als diese ID
      .flatten.node = .frame.id[length(.frame.id)]
      # ...ansonsten...
    } else{
      # ...hole subset aus frame mit letzten (nach dem letzten key liegend) Elementen (! OHNE LeafNode -> ...:(length(.frame) - 1)... )...
      .flatten.node = .frame[(.knm.max + 1):(length(.frame) - 1)]
      # ...falls die Laufnummer von den Zwischenknoten != 0 ist (es gab also parallel noch weitere Knoten mit dieser ID), dann haenge die Laufnummer an (der Wert wird also in den Namen uebernommen)...
      .flatten.node = if_else(.flatten.node == 0, names(.flatten.node), paste0(names(.flatten.node), '_', .flatten.node))
      # ...verbinde die Knotennamen mit '.' & haenge LeafNode-ID mit an
      .flatten.node = paste0(paste0(.flatten.node, collapse = '.'), '.', .frame.id[length(.frame.id)])
    }
      
    # benutze flatten.node als Name des letzen LeafNode-Wertes (.value)
    .value = paste0(.flatten.node, ':', .value)

    # Speichere in LookUp und aggregiere pro Schluessel(-kombi) (keys MIT values) alle 'FlattenNodes' (.value wurden wieder in key:value string gewandelt, so dass hier das unique greift)
    # Hier werden also nach den spezifischen VALUES der Schluessel die flattenNodes einsortiert (bsp.: "Patient_ID:6923 # Tumor_ID:992827 # OP:0")
    Aggregate.LookUp[[.tab.key]] = unique(c(Aggregate.LookUp[[.tab.key]], .value))
  } else
    # ...ansonsten gib Warning aus
    warning('in PATH:\n\'', .x, '\'\n!!! last node ist key node, please check')
})
# Ueberpruefe Ausgabe #
# ls(Aggregate.LookUp)


# Erzeuge LookUp
Aggregate_2.LookUp = new.env(hash = TRUE)

# iteriere ueber Schluessel(-kombis)...
walk(sort(ls(Aggregate.LookUp)), ~{
  # Knoten sind durch ' # ' getrennt und und die Werte haengen mit ':' an den jeweiligen Knoten an (key:value # key:value # key:value # .....)
  
  # zerlege den keystring in key-value Liste
  .key.vec = unlist(strsplit(.x, ' # ')) %>% strsplit(':')
  # erzeuge named vektor mit keys als Namen
  .frame.key = map_chr(.key.vec, ~.x[[2]])
  names(.frame.key) = map_chr(.key.vec, ~.x[[1]])
  # erzeuge ID_Vektor
  .frame.key.id = names(.frame.key)
  
  # zerlege den valuestring in key-value Liste (Knoten liegen im Vektor einzeln vor, KEIN konkatenierter string mit ' # ' Trenner)
  .value.vec = strsplit(Aggregate.LookUp[[.x]], ':')

  # erzeuge named vector mit keys als Namen
  .frame.value = map_chr(.value.vec, ~.x[[2]])
  names(.frame.value) = map_chr(.value.vec, ~.x[[1]])
  .frame.value.id = names(.frame.value)
  
  # erzeuge Schluessel aus den IDs des aktuellen Schluessels (.x)
  .tab.key = paste0(.frame.key.id, collapse = ' # ')
  
  # jetzt werden nach den allgemeinen IDs der Schluessel die Werte einsortiert, diese ergeben dann die spezifischen Tabellen
  # (bsp.: "Patient_ID # Tumor_ID # OP" (NUR keys KEINE values))
  Aggregate_2.LookUp[[.tab.key]] = c(Aggregate_2.LookUp[[.tab.key]], list(c(.frame.key, .frame.value)))
})

# ueberpruefe Ausgabe #
# ls(Aggregate_2.LookUp)
# Aggregate_2.LookUp[['Patient_ID']]
# bind_rows(Aggregate_2.LookUp[['Patient_ID']])

# erstelle LookUp
Table.LookUp = new.env(hash = TRUE)

# iteriere ueber allgemeine Schluessel
walk(sort(ls(Aggregate_2.LookUp)), ~{
  # ...gib aktuellen Schluessel aus
  print(.x)
  # erzeuge Tabellen mit dplyr::bind_rows
  Table.LookUp[[.x]] = dplyr::bind_rows(Aggregate_2.LookUp[[.x]])
})

# Ausgabe der Tabellen-IDs in den Hashes
print(sort(ls(Table.LookUp)))

# gib Tabellen aus
walk(sort(ls(Table.LookUp)), ~{
  print(.x)
  print(Table.LookUp[[.x]])
})






# Rekursion ueber eine Liste von Listen, wie man sie von 'xml2::as_list(xml2::read_xml(<file.xml>))' erhaelt
# Es werden keine Knoten ausgelassen
# Knoten werden durchnummeriert und Pfade werden mit Endknoten-Wert ausgegeben
recListdefault = function(.list, .actual.id = '', .enum.id = '', .add.id.to.path = FALSE, .path = NULL, .level = 1){
  # alle Attributswerte (named vector)
  .atts = attributes(.list)
  # alle Attributsnamen
  .atts.ids = names(.atts)
  
  # Falls es weitere (Sub-)Knoten gibt...
  if('names' %in% .atts.ids){
    # Falls der Knoten in den Pfad eingebracht werden soll (da es weitere Subknoten gibt (-> 'names), wird der nummerierte Knotennamen (.enum.id) eingebracht)
    if(.add.id.to.path)
      .path = c(.path, paste0(.actual.id, ':', .enum.id))
    
    # Attribut 'names' wird entfernt...
    .atts.ids = .atts.ids[!.atts.ids %in% c('names')]
    
    # Falls es weitere Attribute gibt...
    if(length(.atts.ids) > 0){
      # ...baue key:value vector auf
      .key.value = sort(paste0(.atts.ids, ':', .atts[.atts.ids]))
      # ...und pflege in den Pfad ein
      .path = c(.path, .key.value)
    }
    
    # hole alle Knotennamen aus Attribut 'names'
    .ids = names(.list)
    
    # Erzeuge Nummerierungsvektor
    .ids.enum = enumEntsInVec(.ids, .append = FALSE, .single.tag = '0')
    
    # setze "default Knoten darf in Pfad"
    .add.id.to.path = TRUE
    
    # iteriere ueber Knotennamen (Positionen im Knotennamenvektor)
    for(.pos in seq_along(.ids))
      recListdefault(.list[[.pos]], .ids[[.pos]], .ids.enum[[.pos]], .add.id.to.path, .path, .level + 1)
    # Falls wir am Ende eines Knotens sind...
  } else{
    # ...setze die (Value-)Liste zu einem String zusammen (Kommasepariert)
    # ...falls der Endknoten alleine stand
    if(.enum.id %in% '0')
      .key.value = paste0(.actual.id, ':', paste0(sort(unlist(.list)), collapse = ', '))
    # ...falls es den Endknoten oefters gab, haenge die Nummerierung feste am Namen an
    else
      .key.value = paste0(paste0(.actual.id, '_', .enum.id), ':', paste0(sort(unlist(.list)), collapse = ', '))
    
    # ...setze den Pfad-String zusammen (Pfad bis hierher)
    .path.string = paste0(paste0(.path, collapse = ' # '), ' # ', .key.value)
    
    # Falls es diese Eintrag schon einmal gab...
    if(env_has(LookUpDefault, .path.string)){
      # ...und die Werte unterschiedlich sind, breche ab
      if(.value != LookUpDefault[[.path.string]])
        stop(paste0('PATH: ', .path.string, '\nOLD VALUE: ', LookUpDefault[[.path.string]], '\nNEW VALUE: ', .value))
    }
    
    # ...ansonsten speichere (Pfad + aktueller Knotennamen (ohne Nummerierung)) + Wert)
    LookUpDefault[[.path.string]] = 1L
    
    # Falls es weitere Attribute in den LeafNodes gibt (Attribut 'names' war im else-Block nicht dabei)...
    if(length(.atts.ids) > 0){
      # ...baue key:value vector auf, ergaenze keys mit aktueller Knoten_ID (um zugehoerigkeit zu zeigen; -> <KnotenID>_<key>:<value>)
      .key.value = paste0(.actual.id, '_', .atts.ids, ':', .atts[.atts.ids])
      
      walk(.key.value, ~{
        .path.string = paste0(paste0(.path, collapse = ' # '), ' # ', .x)
        LookUpDefault[[.path.string]] = 1L
      })
    }
  }
}

# Erstelle Lookup
LookUpDefault = new.env(hash = TRUE)
# Laufe Rekursion
recListdefault(lxml)
# Gib Pfade aus
print(sort(ls(LookUpDefault)))





# Standard Rekursion ueber eine Liste von Listen, wie man sie von 'xml2::as_list(xml2::read_xml(<file.xml>))' erhaelt
# Es findet keine Veraenderung der Knoten statt, alle Pfade werden ohne Werte ausgegeben (uniqueness) sowie
# alle Endpunkt-Werte, Attribute und Knotennamen
recListGetNodesPathsAttributes = function(.list, .actual.id = '', .add.id.to.path = FALSE, .path = NULL, .level = 1){
  # alle Attributswerte (named vector)
  .atts = attributes(.list)
  # alle Attributsnamen
  .atts.ids = names(.atts)
  
  # Falls es weitere (Sub-)Knoten gibt...
  if('names' %in% .atts.ids){
    # Falls der Knoten in den Pfad eingebracht werden soll
    if(.add.id.to.path)
      .path = c(.path, .actual.id)
    
    # Attribut 'names' wird entfernt...
    .atts.ids = .atts.ids[!.atts.ids %in% c('names')]
    
    # Falls es weitere Attribute gibt...
    if(length(.atts.ids) > 0){
      walk(.atts.ids, ~{
        Attributes.LookUp[[.x]] = 1L
      })
    }

    # hole alle Knotennamen aus Attribut 'names'
    .ids = names(.list)

    # setze "default Knoten darf in Pfad"
    .add.id.to.path = TRUE
    
    # iteriere ueber Knotennamen (Positionen im Knotennamenvektor)
    for(.pos in seq_along(.ids)){
      Nodes.LookUp[[.ids[[.pos]]]] = 1L
      recListGetNodesPathsAttributes(.list[[.pos]], .ids[[.pos]], .add.id.to.path, .path, .level + 1)
    }
    # Falls wir am Ende eines Knotens sind...
  } else{
    # ...setze die (Value-)Liste zu einem String zusammen (Kommasepariert)
    .value = paste0(sort(unlist(.list)), collapse = ', ')
    # ...setze den Pfad-String zusammen (Pfad bis hierher + aktueller Knotenname + Wert
    .path.string = paste0(paste0(.path, ' # ', collapse = ''), .actual.id, ' : <VALUE>')
    
    Paths.LookUp[[.path.string]] = 1L
    Values.LookUp[[.value]] = 1L
    Level.vec <<- c(Level.vec, .level)
    
    # Falls es weitere Attribute gibt...
    if(length(.atts.ids) > 0){
      walk(.atts.ids, ~{
        Attributes.LookUp[[.x]] = 1L
      })
    }
  }
}

# Erzeuge LookUp Tabellen
Paths.LookUp = new.env(hash = TRUE)
Nodes.LookUp = new.env(hash = TRUE)
Values.LookUp = new.env(hash = TRUE)
Attributes.LookUp = new.env(hash = TRUE)
Level.vec = integer(0)

# Laufe Rekursion
recListGetNodesPathsAttributes(lxml)

# Gib LookUp aus
sort(ls(Paths.LookUp))
sort(ls(Nodes.LookUp))
sort(ls(Values.LookUp))
sort(ls(Attributes.LookUp))
table(Level.vec, useNA = 'ifany')


# Ausspielen der Einzeltabellen ####
for (i in names(Table.LookUp)){
  print(i)
  if(i == 'Patient_ID'){
    Tab_Patient_ID <- get(i, envir = Table.LookUp)
  } else{
    position_last_hashtag <- max(unlist(stringr::str_locate_all(i, "#")))
    new_name <- paste0('Tab_', substr(i, position_last_hashtag+2, nchar(i)))
    assign(x = new_name, value =  get(i, envir = Table.LookUp))
    rm(position_last_hashtag, new_name)
  }
}

# nur Tabellen behalten
gdata::keep(list = objects(pattern = '^Tab_'), sure = T)
