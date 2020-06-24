(* ::Package:: *)

<< KerrGeodesics`;
<< ApplicationTools`;

packageDirectory = NotebookDirectory[]; 

packages =
{ 
  "ConstantsOfMotion",
  "KerrGeoOrbit",
  "OrbitalFrequencies",
  "SpecialOrbits",
  "NearHorizonGeoOrbit"
};

packageSymbols = Map[# -> DocumentedSymbols["KerrGeodesics", #] &, packages];

undocumentedSymbols = Map[# -> UndocumentedSymbols["KerrGeodesics", #] &, packages] /. (_ -> {}) -> Sequence[];
Map[Print["Undocumented symbols for package "<>#[[1]]<>" skipped:\n", #[[2]]]&, undocumentedSymbols];

Print["Building symbol reference pages"];
docPackage[package_ -> symbols_] :=
  Map[(Print[#]; BuildSymbolReference["KerrGeodesics", #, "Source", FileNameJoin[{Directory[],"Documentation","English","ReferencePages","Symbols",#}]<>".nb"]) &, symbols];
Scan[docPackage, packageSymbols];

Print["Building guides"];
sourceGuides = FileNames["*.md", FileNameJoin[{"Source", "Documentation", "English", "Guides"}], Infinity];
destGuides =
  FileNameJoin[{Directory[], FileNameDrop[DirectoryName[#], 1],
      FileBaseName[#] <> ".nb"}] & /@ sourceGuides;
MapThread[BuildGuide, {sourceGuides, destGuides}];

Print["Building tutorials"];
tutorialSources = FileNames["*.md", FileNameJoin[{"Source", "Documentation", "English", "Tutorials"}], Infinity];
Map[(Print[#]; BuildTutorial[FileNameJoin[{Directory[], #}], FileNameJoin[{Directory[],"Documentation", "English","Tutorials",FileBaseName[#]}]<>".nb"])&, tutorialSources];

Print["Indexing Documentation"];
BuildIndex["KerrGeodesics"];

Print["Done"];
