unit Unit1;

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, ExtCtrls, Jpeg;

type
  Tpt = record
    // coordonnées World
    x, y, z: real;
    // coordonnées observateur
    xo, yo, zo: real;
    // coordonnées écran
    x2d, y2d: integer;
    // Eloignement moyen (= zO moyen)
    zmoyen: real;
    // Equation du plan de la face (de type Ax+By+Cz+D=0)
    A, B, C, D: real;
    ID: integer;
  end;

  TFace = array of Tpt;
  TObj = array of TFace;

  TForm1 = class(TForm)
    PaintBox1: TPaintBox;
    btn_Affichage: TButton;
    scb_h: TScrollBar;
    scb_v: TScrollBar;
    Label1: TLabel;
    Label2: TLabel;
    btn_Charger: TButton;
    OpenDialog1: TOpenDialog;
    Memo1: TMemo;
    edt_Distance: TEdit;
    edt_Zoom: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    chk_Repere: TCheckBox;
    btn_AfficheFace: TButton;
    edt_FaceNum: TEdit;
    btn_Intersect: TButton;
    edt_intersect1: TEdit;
    edt_intersect2: TEdit;
    chk_Ajuster: TCheckBox;
    btn_FacesShow: TButton;
    btn_MemoClear: TButton;
    Label5: TLabel;
    Label6: TLabel;
    btn_Verifier: TButton;
    btn_TC: TButton;
    chk_Ortho: TCheckBox;
    btn_JPEG: TButton;
    SaveDialog1: TSaveDialog;
    scb_w: TScrollBar;
    Label7: TLabel;
    Label8: TLabel;
    procedure btn_AffichageClick(Sender: TObject);
    procedure DoChangeAngle(Sender: TObject);
    procedure DoCreate(Sender: TObject);
    procedure btn_ChargerClick(Sender: TObject);
    procedure btn_AfficheFaceClick(Sender: TObject);
    procedure btn_IntersectClick(Sender: TObject);
    procedure btn_FacesShowClick(Sender: TObject);
    procedure btn_MemoClearClick(Sender: TObject);
    procedure DoMouseUp(Sender: TObject; Button: TMouseButton;
      Shift: TShiftState; x, y: integer);
    procedure btn_VerifierClick(Sender: TObject);
    procedure btn_TCClick(Sender: TObject);
    procedure DoOrthoChange(Sender: TObject);
    procedure btn_JPEGClick(Sender: TObject);
    procedure DoChangeParam(Sender: TObject);
    procedure DoMouseDown(Sender: TObject; Button: TMouseButton;
      Shift: TShiftState; X, Y: Integer);
    procedure DoMouseMove(Sender: TObject; Shift: TShiftState; X, Y: Integer);
    procedure DoMemoChange(Sender: TObject);
    procedure DoWheelUp(Sender: TObject; Shift: TShiftState; MousePos: TPoint;
      var Handled: Boolean);
    procedure DoWheelDown(Sender: TObject; Shift: TShiftState; MousePos: TPoint;
      var Handled: Boolean);
  private
    { Déclarations privées }
  public
    { Déclarations publiques }

  end;

var
  Form1: TForm1;
  obj: TObj;
  face: TFace;
  som: array of Tpt;
  nbsom: integer;
  nbface: integer;
  sinus: array [0 .. 359] of real;
  cosinus: array [0 .. 359] of real;
  zoom: integer;
  distance: integer;
  TmpBmp: TBitmap;
  calculated : boolean;
  MouseIsDown : boolean;
  PDown : TPoint;
  PActually : TPoint;

implementation

const
  pi = 3.14159;
  tolerance = 0.01;
{$R *.DFM}

procedure debug(s: string);
begin
  Form1.Memo1.Lines.Add(s);
end;

procedure CalculPoint(var p: Tpt; angle_h, angle_v, angle_w, centrex, centrey: integer);
var
  savex : integer;
begin
  with p do
  begin
    // système world --> système de l'observateur
    zo := -x * (cosinus[angle_h] * cosinus[angle_v]) - y *
      (sinus[angle_h] * cosinus[angle_v]) - z * sinus[angle_v] + distance;
    xo := (-x * sinus[angle_h] + y * cosinus[angle_h]);
    yo := (-x * (cosinus[angle_h] * sinus[angle_v]) - y *
      (sinus[angle_h] * sinus[angle_v]) + z * cosinus[angle_v]);

    // projection
    if Form1.chk_Ortho.Checked then
    begin
      x2d := round(xo * distance / 100.0 * zoom / 512.0) + centrex;
      y2d := centrey - (round(yo * distance / 100.0 * zoom / 512.0));
      savex := x2d;
      x2d := round(cosinus[angle_w]*(x2d-centrex) - sinus[angle_w]*(centrey-y2d)) + centrex;
      y2d := centrey - round(sinus[angle_w]*(savex-centrex) + cosinus[angle_w]*(centrey-y2d));
    end
    else
    begin
      x2d := round(xo * zoom / zo) + centrex;
      y2d := centrey - (round(yo * zoom / zo));
      savex := x2d;
      x2d := round(cosinus[angle_w]*(x2d-centrex) - sinus[angle_w]*(centrey-y2d)) + centrex;
      y2d := centrey - round(sinus[angle_w]*(savex-centrex) + cosinus[angle_w]*(centrey-y2d));
      {x' = cos(theta)*(x-xc) - sin(theta)*(y-yc) + xc
      y' = sin(theta)*(x-xc) + cos(theta)*(y-yc) + yc}
    end;

  end;
end;

procedure TriFacesSimple;
// tir des faces en fontion de l'éloignement du centre de la face
// par rapport à l'observateur
// on travaille dans le système de l'observateur
var
  curface, cursommet: integer;
  tmpface: TFace;
  zmoyen: real;
  permut: boolean;
begin
  // calcule le zmoyen (= éloignement) pour chaque face
  for curface := 0 to length(obj) - 1 do
  begin
    zmoyen := 0;
    for cursommet := 0 to length(obj[curface]) - 1 do
    begin
      zmoyen := zmoyen + obj[curface][cursommet].zo;
    end;
    // éloignement stocké dans 1er point de la face
    obj[curface][0].zmoyen := zmoyen / length(obj[curface]);
  end;

  // tri à bulle
  repeat
    permut := false;
    for curface := 0 to length(obj) - 2 do
    begin
      // si face n+1 est plus éloignée que face n ==> permutation
      if obj[curface + 1][0].zmoyen > obj[curface][0].zmoyen then
      begin
        permut := true;
        tmpface := obj[curface];
        obj[curface] := obj[curface + 1];
        obj[curface + 1] := tmpface;
      end;
    end;
  until permut = false;
  // fin tri à bulle.
end;

function Signe(r: real): integer;
begin
  if r >= 0 then
    Signe := 1
  else
    Signe := -1;
end;

function Arrondi(r: real): real;
begin
  if abs(r) < tolerance then
    Arrondi := 0
  else
    Arrondi := r;
end;

function LimiteXY(p, Q: TFace): boolean;
// test les limites de projections 2D
// ==> si les boites englobantes ne se chevauchent pas --> renvoie true
var
  maxPx, maxPy, minPx, minPy, maxQx, maxQy, minQx, minQy, i: integer;
begin
  // calule min et max x et y de P
  maxPx := p[0].x2d;
  minPx := p[0].x2d;
  maxPy := p[0].y2d;
  minPy := p[0].y2d;
  for i := 1 to length(p) - 1 do
  begin
    if p[i].x2d > maxPx then
      maxPx := p[i].x2d;
    if p[i].x2d < minPx then
      minPx := p[i].x2d;
    if p[i].y2d > maxPy then
      maxPy := p[i].y2d;
    if p[i].y2d < minPy then
      minPy := p[i].y2d;
  end;
  // calule min et max x et y de Q
  maxQx := Q[0].x2d;
  minQx := maxQx;
  maxQy := Q[0].y2d;
  minQy := maxQy;
  for i := 1 to length(Q) - 1 do
  begin
    if Q[i].x2d > maxQx then
      maxQx := Q[i].x2d;
    if Q[i].x2d < minQx then
      minQx := Q[i].x2d;
    if Q[i].y2d > maxQy then
      maxQy := Q[i].y2d;
    if Q[i].y2d < minQy then
      minQy := Q[i].y2d;
  end;
  // teste le chevauchement
  if (maxPx < minQx) or (maxQx < minPx) or (maxPy < minQy) or
    (maxQy < minPy) then
    LimiteXY := true // disjoints
  else
    LimiteXY := false; // chevauchement
end;

function Devant(p, Q: TFace; angle_h, angle_v: integer): integer;
// tous les point de Q sont-ils devant le plan P ?
var
  Ap, Bp, Cp, Dp: real;
  i, positif, negatif, cotepoint: integer;
  position, tempo, xobs, yobs, zobs: real;
begin
  Ap := p[0].A;
  Bp := p[0].B;
  Cp := p[0].C;
  Dp := p[0].D;

  // positon observateur dans le système world;
  xobs := Arrondi(distance * cosinus[angle_h] * cosinus[angle_v]);
  yobs := Arrondi(distance * sinus[angle_h] * cosinus[angle_v]);
  zobs := Arrondi(distance * sinus[angle_v]);

  // De quel coté du plan P sont les points de la face Q ?
  positif := 0;
  negatif := 0;
  cotepoint := 0;
  for i := 0 to length(Q) - 1 do
  begin
    tempo := Arrondi(Ap * Q[i].x + Bp * Q[i].y + Cp * Q[i].z + Dp);
    if tempo > 0 then
      inc(positif);
    if tempo < 0 then
      inc(negatif);
    // on ignore le cas ou tempo = 0
    // ==> ce point de Q est dans plan P
    // ==> neutralité
  end;

  // points de Q de chaque coté de P
  if (positif > 0) and (negatif > 0) then
  begin
    // on ne peut rien déterminer. Q "à cheval" sur le plan Pbegin
    Devant := 0;
    // debug('fonction Devant : P à cheval sur le plan Q');
    Exit;
  end
  // tous les point de Q sont du même coté du plan P
  else
  begin
    if positif > 0 then
      cotepoint := 1;
    if negatif > 0 then
      cotepoint := -1;
    // dans quel demi espace est l'observateur ?
    position := Signe(Arrondi(Ap * xobs + Bp * yobs + Cp * zobs + Dp));
    if position = cotepoint then // Q du même coté que l'observateur ?
    begin
      Devant := 1; // oui
      Exit;
    end
    else // Q est derrière P ==> permutation
    begin
      Devant := -1;
      Exit;
    end;
  end;
end;

function Derriere(p, Q: TFace; angle_h, angle_v: integer): integer;
// tous les point de P sont-ils derrière le plan Q ?
var
  Aq, Bq, Cq, Dq: real;
  i, positif, negatif, cotepoint: integer;
  position, tempo, xobs, yobs, zobs: real;
begin
  Aq := Q[0].A;
  Bq := Q[0].B;
  Cq := Q[0].C;
  Dq := Q[0].D;

  // positon observateur dans le système world;
  xobs := Arrondi(distance * cosinus[angle_h] * cosinus[angle_v]);
  yobs := Arrondi(distance * sinus[angle_h] * cosinus[angle_v]);
  zobs := Arrondi(distance * sinus[angle_v]);

  // De quel coté du plan Q sont les points de la face P ?
  positif := 0;
  negatif := 0;
  cotepoint := 0;
  for i := 0 to length(p) - 1 do
  begin
    tempo := Arrondi(Aq * p[i].x + Bq * p[i].y + Cq * p[i].z + Dq);
    if tempo > 0 then
      inc(positif);
    if tempo < 0 then
      inc(negatif);
    // on ignore le cas ou tempo = 0
    // ==> ce point de Q est dans plan P
    // ==> neutralité
  end;

  // points de Q de chaque coté de P
  if (positif > 0) and (negatif > 0) then
  begin
    // on ne peut rien déterminer. Q "à cheval" sur le plan P
    Derriere := 0;
    // debug('fonction Derriere : Q à cheval sur le plan P');
    Exit;
  end
  else
  // tous les point de P sont du même coté du plan Q
  begin
    if positif > 0 then
      cotepoint := 1;
    if negatif > 0 then
      cotepoint := -1;
    // dans quel demi espace est l'observateur ?
    position := Signe(Arrondi(Aq * xobs + Bq * yobs + Cq * zobs + Dq));
    // les points de P sont du même coté que l'observateur ?
    if position <> cotepoint then
      // oui
      Derriere := 1
    else
      // non permutation
      Derriere := -1;
  end;
end;

function LimiteZ(p, Q: TFace): integer;
var
  i: integer;
  minZP, maxZP, maxZQ, minZQ: real;
begin
  minZP := p[0].zo;
  maxZP := p[0].zo;
  for i := 1 to length(p) - 1 do
  begin
    if p[i].zo < minZP then
      minZP := p[i].zo;
    if p[i].zo > maxZP then
      maxZP := p[i].zo;
  end;
  minZQ := Q[0].zo;
  maxZQ := Q[0].zo;
  for i := 1 to length(Q) - 1 do
  begin
    if Q[i].zo < minZQ then
      minZQ := Q[i].zo;
    if Q[i].zo > maxZQ then
      maxZQ := Q[i].zo;
  end;
  // teste les maxi-mini
  if minZP > maxZQ then
    LimiteZ := 1 // P plus loin que Q
  else if minZQ > maxZP then
    LimiteZ := -1 // Q plus loin que P
  else
    LimiteZ := 0; // indéterminé : chevauchement en Z
end;

function LimiteRegion(p, Q: TFace): boolean;
var
  rP, rQ, r: HRGN;
  ptabP, ptabQ: array of TPoint;
  i: integer;
begin
  // préparation des tableaux de points
  SetLength(ptabP, length(p));
  for i := 0 to length(p) - 1 do
    ptabP[i] := point(p[i].x2d, p[i].y2d);
  SetLength(ptabQ, length(Q));
  for i := 0 to length(Q) - 1 do
    ptabQ[i] := point(Q[i].x2d, Q[i].y2d);
  // creation des régions
  rP := CreatePolygonRgn(ptabP[0], length(p), WINDING);
  rQ := CreatePolygonRgn(ptabQ[0], length(Q), WINDING);
  r := CreatePolygonRgn(ptabQ[0], length(Q), WINDING); // pour init. R
  // intersection
  i := CombineRgn(r, rP, rQ, RGN_AND);
  { case i of
    NULLREGION : debug('combinaison = null region') ;
    ERROR : debug('combinaison = error') ;
    SIMPLEREGION : debug('combinaison = simple') ;
    COMPLEXREGION : debug('combinaison = complexe') ;
    end; }
  if i = NULLREGION then
    LimiteRegion := true
  else
    LimiteRegion := false;
  DeleteObject(rP);
  DeleteObject(rQ);
  DeleteObject(r);
end;

function Eloignement(p, Q: TFace): boolean;
// mesure l'éloignemetLimiteRegion du centre des 2 faces (= moyenne des zo)
// renvoie true si moyenne de P > moyenne de Q
var
  zmoyenP, zmoyenQ: real;
  i: integer;
begin
  zmoyenP := 0;
  for i := 0 to length(p) - 1 do
    zmoyenP := zmoyenP + p[i].zo;
  zmoyenP := zmoyenP / length(p);

  zmoyenQ := 0;
  for i := 0 to length(Q) - 1 do
    zmoyenQ := zmoyenQ + Q[i].zo;
  zmoyenQ := zmoyenQ / length(Q);

  if zmoyenP > zmoyenQ then
    Eloignement := true
  else
    Eloignement := false;
end;

procedure Ajuster;
var
  i, j: integer;
  xmin, xmax, ymin, ymax, zmin, zmax, xdelta, ydelta, zdelta: real;
  ratio, K: real;
begin
  K := 120;
  xmin := obj[0][0].x;
  xmax := obj[0][0].x;
  ymin := obj[0][0].y;
  ymax := obj[0][0].y;
  zmin := obj[0][0].z;
  zmax := obj[0][0].z;

  for i := 0 to length(obj) - 1 do
    for j := 0 to length(obj[i]) - 1 do
    begin
      if xmin < obj[i][j].x then
        xmin := obj[i][j].x;
      if xmax > obj[i][j].x then
        xmax := obj[i][j].x;
      if ymin < obj[i][j].y then
        ymin := obj[i][j].y;
      if ymax > obj[i][j].y then
        ymax := obj[i][j].y;
      if zmin < obj[i][j].z then
        zmin := obj[i][j].z;
      if zmax > obj[i][j].z then
        zmax := obj[i][j].z;
    end;
  xdelta := (xmax + xmin) / 2;
  ydelta := (ymax + ymin) / 2;
  zdelta := (zmax + zmin) / 2;

  if xmax - xmin > ymax - ymin then
    ratio := xmax - xmin
  else
    ratio := ymax - ymin;
  if zmax - zmin > ratio then
    ratio := zmax - zmin;
  ratio := abs(K / ratio);

  for i := 0 to length(obj) - 1 do
    for j := 0 to length(obj[i]) - 1 do
    begin
      obj[i][j].x := ratio * (obj[i][j].x - xdelta);
      obj[i][j].y := ratio * (obj[i][j].y - ydelta);
      obj[i][j].z := ratio * (obj[i][j].z - zdelta);
    end;
end;

procedure AfficheFace(face: integer);
var
  i: integer;
  f: array of TPoint;
begin
  SetLength(f, length(obj[face]));
  for i := 0 to length(obj[face]) - 1 do
  begin
    f[i].x := obj[face][i].x2d;
    f[i].y := obj[face][i].y2d;
  end;
  with Form1.PaintBox1.Canvas do
  begin
    Brush.Color := clGreen;
    Pen.Color := clLtGray;
    Polygon(f);
  end;
end;

function ZBufferIntersect(p, Q: TFace; centrex, centrey: integer): boolean;
// Détermine  la face (P ou Q) la plus éloignée par le Zbuffer
// ==> compare les Z (système observateur) du centtre
// de la région intersection des régions de P et Q
var
  r, rP, rQ: HRGN;
  rect: TRect;
  ptabP, ptabQ: Array of TPoint;
  maxPx, maxPy, minPx, minPy, maxQx, maxQy, minQx, minQy, i, j: integer;
  zP, zQ: real;
  centre: TPoint;
  x1, y1, z1, x2, y2, z2, x3, y3, z3: real;
  Pa, Pb, Pc, Pd, Qa, Qb, Qc, Qd: real;
  nbx, nby, nbin, nb: integer;

begin
  // si les région sont disjointe pas d'intersection
  if LimiteRegion(p, Q) then
  begin
    ZBufferIntersect := true;
    Exit;
  end;
  // Création des régions (projection 2D de P et Q)
  SetLength(ptabP, length(p));
  for i := 0 to length(p) - 1 do
    ptabP[i] := point(p[i].x2d, p[i].y2d);
  SetLength(ptabQ, length(Q));
  for i := 0 to length(Q) - 1 do
    ptabQ[i] := point(Q[i].x2d, Q[i].y2d);
  rP := CreatePolygonRgn(ptabP[0], length(p), WINDING);
  rQ := CreatePolygonRgn(ptabQ[0], length(Q), WINDING);
  r := CreatePolygonRgn(ptabQ[0], length(Q), WINDING); // pour init. R
  // Création de l'intersection
  i := CombineRgn(r, rP, rQ, RGN_AND);
  { case i of
    NULLREGION:
    debug('combinaison = null region');
    ERROR:
    debug('combinaison = error');
    SIMPLEREGION:
    debug('combinaison = simple');
    COMPLEXREGION:
    debug('combinaison = complexe');
    end; }

  if i = NULLREGION then
    debug('Null Region');
  if i = ERROR then
    debug('ERROR Region');

  // Affichage de la région intersection de P et Q
  FrameRgn(Form1.PaintBox1.Canvas.Handle, r,
    Form1.PaintBox1.Canvas.Brush.Handle, 2, 2);

  // r = région intersection.
  // calcul du centre de la région r
  GetRgnBox(r, rect);
  nbx := 0;
  nby := 0;
  nbin := 0;
  nb := 0;
  for i := rect.left to rect.right do
    for j := rect.top to rect.Bottom do
    begin
      inc(nb);
      if PtInRegion(r, i, j) then
      begin
        inc(nbin);
        nbx := nbx + i;
        nby := nby + j;
      end;
    end;

  DeleteObject(r);
  DeleteObject(rP);
  DeleteObject(rQ);

  if nbin = 0 then
  begin
    ZBufferIntersect := true;
    Exit;
  end;

  centre.x := nbx div nbin;
  centre.y := nby div nbin;

  // affichage du centre de l'intersection
  rect.left := centre.x - 2;
  rect.right := centre.x + 2;
  rect.top := centre.y - 2;
  rect.Bottom := centre.y + 2;
  Form1.PaintBox1.Canvas.FrameRect(rect);

  x1 := p[0].xo;
  y1 := p[0].yo;
  z1 := p[0].zo;
  x2 := p[1].xo;
  y2 := p[1].yo;
  z2 := p[1].zo;
  x3 := p[2].xo;
  y3 := p[2].yo;
  z3 := p[2].zo;
  Pa := Arrondi(y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2));
  Pb := Arrondi(-x1 * (z2 - z3) + x2 * (z1 - z3) - x3 * (z1 - z2));
  Pc := Arrondi(x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2));
  Pd := Arrondi(-x1 * (y2 * z3 - y3 * z2) + x2 * (y1 * z3 - y3 * z1) - x3 *
    (y1 * z2 - y2 * z1));
  // debug(' de P0 par  calcul = ' + FloatToStr(-(Pa * p[0].xo + Pb * p[0].yo +
  // Pd) / Pc));
  // debug(' de P0 initial = ' + FloatToStr(z1));
  x1 := Q[0].xo;
  y1 := Q[0].yo;
  z1 := Q[0].zo;
  x2 := Q[1].xo;
  y2 := Q[1].yo;
  z2 := Q[1].zo;
  x3 := Q[2].xo;
  y3 := Q[2].yo;
  z3 := Q[2].zo;
  Qa := Arrondi(y1 * (z2 - z3) + y2 * (z3 - z1) + y3 * (z1 - z2));
  Qb := Arrondi(-x1 * (z2 - z3) + x2 * (z1 - z3) - x3 * (z1 - z2));
  Qc := Arrondi(x1 * (y2 - y3) - x2 * (y1 - y3) + x3 * (y1 - y2));
  Qd := Arrondi(-x1 * (y2 * z3 - y3 * z2) + x2 * (y1 * z3 - y3 * z1) - x3 *
    (y1 * z2 - y2 * z1));

  zP := -Pd / ((Pa * (centre.x - centrex) / zoom) + (Pb * (centrey - centre.y) /
    zoom) + Pc);
  zQ := -Qd / ((Qa * (centre.x - centrex) / zoom) + (Qb * (centrey - centre.y) /
    zoom) + Qc);

  if zP >= zQ then
    ZBufferIntersect := true
  else
    ZBufferIntersect := false;
end;

procedure TestComplet(angle_h, angle_v, angle_w, centrex, centrey: integer);
// tri des face suivant leurposition dans l'espace
// par raport à l'observateur
// on travaille dans le système World.
type
  binome = array [1 .. 2] of integer;
var
  curface, i, j, nbcycle, nbpermut, nbpermutb, savnbpermutb,
    nbiteration, nbindetermine: integer;
  tmpface: TFace;
  x1, y1, z1, x2, y2, z2, x3, y3, z3: real;
  permut: boolean;
  p, Q: TFace;
  tabpermut, tabpermutb: array of binome;

  procedure DoPermut(prio: integer);
  var
    A, B, l: integer;
    trouve: boolean;
  begin
    // détermine les ID des faces P et Q
    A := p[0].ID;
    B := Q[0].ID;
    // A t on déjà faiit cette permutation ?
    trouve := false;
    for l := 0 to nbpermut - 1 do
      if ((tabpermut[l][1] = A) and (tabpermut[l][2] = B)) or
        ((tabpermut[l][1] = B) and (tabpermut[l][2] = A)) then
        trouve := true;
    if (trouve) and (prio = 0) then
    begin
      // debug('permut bloqué !!');
      inc(nbpermutb);
      SetLength(tabpermutb, nbpermutb);
      tabpermutb[nbpermutb - 1][1] := p[0].ID;
      tabpermutb[nbpermutb - 1][2] := Q[0].ID;
      Exit; // oui : sortie (pas de permutation)
    end
    else
    // non, on permute
    begin
      // mise à jour du tableau
      inc(nbpermut);
      SetLength(tabpermut, nbpermut);
      if p[0].ID < Q[0].ID then
      begin
        tabpermut[nbpermut - 1][1] := p[0].ID;
        tabpermut[nbpermut - 1][2] := Q[0].ID;
      end
      else
      begin
        tabpermut[nbpermut - 1][1] := Q[0].ID;
        tabpermut[nbpermut - 1][2] := p[0].ID;
      end;
      // debug(IntToStr(P[0].ID)+ ' '+IntToStr(Q[0].ID));
      // permutation
      permut := true;
      inc(nbpermut);
      tmpface := obj[i];
      obj[i] := obj[j];
      obj[j] := tmpface;
    end;
  end;

begin
  // calcul de l'équation de plan pour chaque face
  // les données sont stockées dans le 1er point de chaque face (obj[curface][0])
  for curface := 0 to length(obj) - 1 do
  begin
    x1 := obj[curface][0].x;
    y1 := obj[curface][0].y;
    z1 := obj[curface][0].z;
    x2 := obj[curface][1].x;
    y2 := obj[curface][1].y;
    z2 := obj[curface][1].z;
    x3 := obj[curface][2].x;
    y3 := obj[curface][2].y;
    z3 := obj[curface][2].z;
    obj[curface][0].A := Arrondi(y1 * (z2 - z3) + y2 * (z3 - z1) + y3 *
      (z1 - z2));
    obj[curface][0].B := Arrondi(-x1 * (z2 - z3) + x2 * (z1 - z3) - x3 *
      (z1 - z2));
    obj[curface][0].C := Arrondi(x1 * (y2 - y3) - x2 * (y1 - y3) + x3 *
      (y1 - y2));
    obj[curface][0].D := Arrondi(-x1 * (y2 * z3 - y3 * z2) + x2 *
      (y1 * z3 - y3 * z1) - x3 * (y1 * z2 - y2 * z1));
  end;
  // tri
  nbcycle := 0;
  SetLength(tabpermut, 0);
  SetLength(tabpermutb, 0);
  nbpermut := 0;
  permut := false;
  nbpermut := 0;
  TriFacesSimple;

  nbpermutb := 0;
  savnbpermutb := 0;
  nbiteration := 0;
  repeat
    inc(nbiteration);
    savnbpermutb := nbpermutb;
    nbpermutb := 0;
    SetLength(tabpermutb, 0);
    for i := 0 to length(obj) - 2 do
    begin

      for j := i + 1 to length(obj) - 1 do
      begin

        inc(nbcycle);
        p := obj[i];
        Q := obj[j];

        // test les limites de projections 2D
        if LimiteXY(p, Q) then
          Continue; // projections disjointes

        // teste les limites de  profondeurs (système observateur)
        // inutile si  TriFacesSimple
        { if LimiteZ(p, Q) = 1 then
          Continue // P "derrière" Q
          else if LimiteZ(p, Q) = -1 then // Q "derrière" P
          begin
          DoPermut(0);
          Continue;
          debug('limiteZ');
          end; }

        // Test des régions
        if LimiteRegion(p, Q) then
          Continue;

        // Tous les points de Q devant P ?
        if Devant(p, Q, angle_h, angle_v) = 1 then
          Continue;

        // Tous les points de P derrière le plan Q ?
        if Derriere(p, Q, angle_h, angle_v) = 1 then
          Continue;

        // Tous les points de P devant le plan Q ?
        if (Devant(Q, p, angle_h, angle_v) = 1) then
        begin
          DoPermut(0);
          Continue;
        end;

        // Tous les points de Q derrière le plan P ?
        if (Derriere(Q, p, angle_h, angle_v) = 1) then
        begin
          DoPermut(0);
          Continue;
        end;

       { if ZBufferIntersect(p, Q, centrex, centrey) = false then
        begin
          DoPermut(1); // prioritaire
          Continue;
        end
        else
          Continue;  }

        //debug('inderterminé');
        inc(nbindetermine);
      end; // for i
    end; // for j
    // debug('nbpermutb = ' + IntToStr(nbpermutb));
    // debug('sav nbpermutb = ' + IntToStr(savnbpermutb));
  until (nbpermutb <= savnbpermutb) and ((nbiteration > 1) or (nbpermutb = 0)) ;
  debug('---------------------');
  debug('Nb. permutation = ' + IntToStr(nbpermut));
  debug('Nb. permutation bloc = ' + IntToStr(nbpermutb));
  debug('Nb. itérations = ' + IntToStr(nbiteration));
  debug('Nb. indeterminés = ' + IntToStr(nbindetermine));
end;

// Affichage systeme de Référence
procedure Reference(angle_h, angle_v, angle_w, centrex, centrey: integer);
var
  p: array of Tpt;
  i: integer;
begin
  SetLength(p, 4);
  with p[0] do
  begin
    x := 0;
    y := 0;
    z := 0;
  end;
  with p[1] do
  begin
    x := 100;
    y := 0;
    z := 0;
  end;
  with p[2] do
  begin
    x := 0;
    y := 100;
    z := 0;
  end;
  with p[3] do
  begin
    x := 0;
    y := 0;
    z := 100;
  end;
  // calcul du système de référence
  for i := 0 to 3 do // pour chaque points
    CalculPoint(p[i], angle_h, angle_v, angle_w, centrex, centrey);

  // affichage système de référence
  with TmpBmp.Canvas do
  begin
    Pen.Width := 2;
    // x
    Pen.Color := clRed;
    Moveto(p[0].x2d, p[0].y2d);
    LineTo(p[1].x2d, p[1].y2d);
    // y
    Pen.Color := clGreen;
    Moveto(p[0].x2d, p[0].y2d);
    LineTo(p[2].x2d, p[2].y2d);
    // z
    Pen.Color := clBlue;
    Moveto(p[0].x2d, p[0].y2d);
    LineTo(p[3].x2d, p[3].y2d);
  end;
end;

procedure TForm1.btn_AffichageClick(Sender: TObject);
var
  i, j: integer;
  angle_h, angle_v, angle_w: integer;
  f: array of TPoint;
  centrex, centrey: integer;
  r: TRect;

begin
  // init.
  TmpBmp := TBitmap.Create;
  TmpBmp.Width := PaintBox1.Width;
  TmpBmp.Height := PaintBox1.Height;
  TmpBmp.Canvas.Lock;

  distance := StrToInt(edt_Distance.Text);
  zoom := StrToInt(edt_Zoom.Text);
  centrex := PaintBox1.Width div 2;
  centrey := PaintBox1.Height div 2;
  angle_h := scb_h.position;
  if angle_h < 0 then
    angle_h := 360 + angle_h;
  angle_v := scb_v.position;
  if angle_v < 0 then
    angle_v := 360 + angle_v;

  angle_w := scb_w.Position;
  if angle_w < 0 then
    angle_w := 360 + angle_w;

  // fond noir
  r.left := 0;
  r.top := 0;
  r.right := PaintBox1.Width;
  r.Bottom := PaintBox1.Height;
  with TmpBmp.Canvas do
  begin
    Brush.Color := clBlack;
    FillRect(r);
  end;

  // Tri des faces !
  if (length(obj)>1) then
  begin
    for i := 0 to length(obj) - 1 do
      for j := 0 to length(obj[i]) - 1 do
        CalculPoint(obj[i][j], angle_h, angle_v, angle_w, centrex, centrey);
  end;

    if not(calculated) and (length(obj)>1) then
    begin
      TestComplet(angle_h, angle_v, angle_w, centrex, centrey);
    end;

  // Affiche les FACES !!!
  TmpBmp.Canvas.Pen.Width := 2;
  for i := 0 to length(obj) - 1 do
  begin
    SetLength(f, length(obj[i]));
    for j := 0 to length(obj[i]) - 1 do
    begin
      f[j].x := obj[i][j].x2d;
      f[j].y := obj[i][j].y2d;
    end;
    with TmpBmp.Canvas do
    begin
      Brush.Color := clWhite;
      Pen.Color := clLtGray;
      Polygon(f);
    end;
  end;
  // Affichage du système de Référence
  if chk_Repere.Checked then
    Reference(angle_h, angle_v, angle_w, centrex, centrey);
  // affichage de la bitmap dans laquelleon a dessiné
  Form1.PaintBox1.Canvas.Draw(0, 0, TmpBmp);
  TmpBmp.Canvas.Unlock;
  TmpBmp.Destroy;
end;

procedure TForm1.btn_ChargerClick(Sender: TObject);
var
  filin: textfile;
  s: string;
  sl: TStringList;
  i, nbpermut: integer;
begin
  nbsom := 0;
  nbface := 0;
  if OpenDialog1.Execute then
  begin
    assignFile(filin, OpenDialog1.FileName);
    reset(filin);
  end
  else
    Exit;
  // entete = sommets
  // lit jusqu'à ligne "sommets" (normalement c'est la 1ere ligne
  // mais pourqoi pas avoir un commentaire avant
  repeat
    readln(filin, s);
    s := LowerCase(trim(s));
  until (pos('sommets', s) = 1) or (eof(filin));
  // si on atteint la fin du fichier sas lire "sommets" ==> erreur
  if eof(filin) then
  begin
    CloseFile(filin);
    Application.MessageBox('Fichier incorrect (pas de section "sommets")',
      'ATTENTION', MB_ICONEXCLAMATION + MB_OK);
    Exit;
  end;
  // lecture des sommets
  repeat
    readln(filin, s);
    s := LowerCase(trim(s));
    if length(s) > 0 then
    begin
      if (pos('faces', s) <> 1) then
      begin
        sl := TStringList.Create();
        sl.Text := StringReplace(s, '.', ',', [rfReplaceAll]);

        sl.Text := StringReplace(s, ' ', #13#10, [rfReplaceAll]);
        inc(nbsom);
        SetLength(som, nbsom);
        som[nbsom - 1].x := StrTOFloat(sl[0]);
        som[nbsom - 1].y := StrTOFloat(sl[1]);
        som[nbsom - 1].z := StrTOFloat(sl[2]);
        if sl.Count <> 3 then
        begin
          Application.MessageBox('Fichier incorrect (nb. de coordonnées <> 3)',
            'ATTENTION', MB_ICONEXCLAMATION + MB_OK);
          sl.Free();
          CloseFile(filin);
          Exit;
        end;
        sl.Free();
      end;
    end;
  until (pos('faces', s) = 1) or (eof(filin));

  // test fin de fichier ==> erreur
  if eof(filin) then
  begin
    CloseFile(filin);
    Application.MessageBox('Fichier incorrect (pas de section "faces")',
      'ATTENTION', MB_ICONEXCLAMATION + MB_OK);
    Exit;
  end;
  // lecture des faces
  repeat
    readln(filin, s);
    if length(s) > 0 then
    begin
      inc(nbface);
      SetLength(obj, nbface);
      sl := TStringList.Create();
      sl.Text := StringReplace(s, ' ', #13#10, [rfReplaceAll]);
      SetLength(obj[nbface - 1], sl.Count);
      if sl.Count > 2 then
      begin
        for i := 0 to (sl.Count - 1) do
        begin
          obj[nbface - 1][i].x := som[StrToInt(sl[i])].x;
          obj[nbface - 1][i].y := som[StrToInt(sl[i])].y;
          obj[nbface - 1][i].z := som[StrToInt(sl[i])].z;
        end; // for
        obj[nbface - 1][0].ID := nbface - 1;
      end // if sl.count ...
      else
      begin
        CloseFile(filin);
        sl.Free;
        Application.MessageBox('Fichier incorrect (nb de sommet incorrect)',
          'ATTENTION', MB_ICONEXCLAMATION + MB_OK);
        Exit;
      end;
      sl.Free;
    end;
    // length > 0
  until eof(filin);
  CloseFile(filin);
  if chk_Ajuster.Checked then
    Ajuster;
  PaintBox1.Invalidate;
end;

procedure TForm1.btn_FacesShowClick(Sender: TObject);
var
  i: integer;
  temps, temps2  : int64;
  r: TRect;
  Color: TColor;
begin
  // fond noir
  r.left := 0;
  r.top := 0;
  r.right := Form1.PaintBox1.Width;
  r.Bottom := Form1.PaintBox1.Height;
  Color := Form1.PaintBox1.Canvas.Brush.Color;
  Form1.PaintBox1.Canvas.Brush.Color := clBlack;
  Form1.PaintBox1.Canvas.FillRect(r);
  Form1.PaintBox1.Canvas.Brush.Color := Color;

  for i := 0 to length(obj) - 1 do
  begin
    QueryPerformanceCounter(temps);
    AfficheFace(i);
    //debug(IntToStr(i));
    repeat
      QueryPerformanceCounter(temps2)
    until (temps2 - temps > 12000);
  end;
end;

procedure TForm1.btn_AfficheFaceClick(Sender: TObject);
begin
  if (edt_FaceNum.Text <> '') and
    (StrToInt(edt_FaceNum.Text) < length(obj)) then
  begin
    AfficheFace(StrToInt(edt_FaceNum.Text));
  end;
end;

procedure TForm1.btn_IntersectClick(Sender: TObject);
var
  r, rP, rQ: TRect;
  maxPx, maxPy, minPx, minPy, maxQx, maxQy, minQx, minQy, i: integer;
  // zP, zQ : real;
  angle_h, angle_v: integer;
  p, Q: TFace;
  tmpbr: Tbrush;
begin
  if length(obj) < 2 then
    Exit;
  if (StrToInt(edt_intersect1.Text) > length(obj) - 1) or
    (StrToInt(edt_intersect2.Text) > length(obj) - 1) then
    Exit;

  p := obj[StrToInt(edt_intersect1.Text)];
  Q := obj[StrToInt(edt_intersect2.Text)];
  // calule min et max x et y de P
  maxPx := p[0].x2d;
  minPx := maxPx;
  maxPy := p[0].y2d;
  minPy := maxPy;
  for i := 1 to length(p) - 1 do
  begin
    if p[i].x2d > maxPx then
      maxPx := p[i].x2d;
    if p[i].x2d < minPx then
      minPx := p[i].x2d;
    if p[i].y2d > maxPy then
      maxPy := p[i].y2d;
    if p[i].y2d < minPy then
      minPy := p[i].y2d;
  end;
  // calule min et max x et y de Q
  maxQx := Q[0].x2d;
  minQx := maxQx;
  maxQy := Q[0].y2d;
  minQy := maxQy;
  for i := 1 to length(Q) - 1 do
  begin
    if Q[i].x2d > maxQx then
      maxQx := Q[i].x2d;
    if Q[i].x2d < minQx then
      minQx := Q[i].x2d;
    if Q[i].y2d > maxQy then
      maxQy := Q[i].y2d;
    if Q[i].y2d < minQy then
      minQy := Q[i].y2d;
  end;
  // interesection des boites limites de P et Q
  rP.left := minPx;
  rP.right := maxPx;
  rP.top := minPy;
  rP.Bottom := maxPy;
  rQ.left := minQx;
  rQ.right := maxQx;
  rQ.top := minQy;
  rQ.Bottom := maxQy;
  IntersectRect(r, rP, rQ);

  angle_h := scb_h.position;
  angle_v := scb_v.position;
  if angle_h < 0 then
    angle_h := 360 + angle_h;
  angle_v := scb_v.position;
  if angle_v < 0 then
    angle_v := 360 + angle_v;

  tmpbr := PaintBox1.Canvas.Brush;
  AfficheFace(StrToInt(edt_intersect1.Text));
  AfficheFace(StrToInt(edt_intersect2.Text));
  debug('Devant(p,q) = ' + IntToStr(Devant(p, Q, angle_h, angle_v)));
  debug('Devant(q,p) = ' + IntToStr(Devant(Q, p, angle_h, angle_v)));
  debug('Derriere(p,q) = ' + IntToStr(Derriere(p, Q, angle_h, angle_v)));
  debug('Derriere(q,p) = ' + IntToStr(Derriere(Q, p, angle_h, angle_v)));

  PaintBox1.Canvas.Brush.Color := clRed;
  PaintBox1.Canvas.FrameRect(r);

  LimiteRegion(p, Q);
  ZBufferIntersect(p, Q, PaintBox1.Width div 2, PaintBox1.Height div 2);

  PaintBox1.Canvas.Brush := tmpbr;
end;

procedure TForm1.btn_JPEGClick(Sender: TObject);
var TmpBmp: TBitmap;
procedure SaveBmpToJpegFile(const Bmp: TBitmap; const FileName: TFileName);
const
  TauxCompressionJpg = 100;
begin
  with TJPEGImage.Create do
  try
    CompressionQuality := TauxCompressionJpg;
    Assign(Bmp);
    SaveToFile(FileName);
  finally
    Free;
  end;
end;
begin
  if Length(obj) = 0 then Exit;

  TmpBmp := TBitmap.Create;
  with TmpBmp do
  try
    Width := PaintBox1.Width;
    Height:= PaintBox1.Height;
    BitBlt(Canvas.Handle, 0, 0, Width, Height, Form1.Canvas.Handle, PaintBox1.Left, PaintBox1.Top, srcCopy);
    if SaveDialog1.Execute then SaveBmpToJpegFile(TmpBmp,saveDialog1.FileName);  // <- avant Free
  finally
    Free;
  end;
end;

procedure TForm1.btn_MemoClearClick(Sender: TObject);
begin
  Memo1.Clear;
end;

procedure TForm1.btn_TCClick(Sender: TObject);
var
  angle_h, angle_v, angle_w, i, j, centrex, centrey: integer;
begin
  angle_h := scb_h.position;
  if angle_h < 0 then
    angle_h := 360 + angle_h;
  angle_v := scb_v.position;
  if angle_v < 0 then
    angle_v := 360 + angle_v;
  angle_w := scb_w.Position;
  if angle_w < 0 then
  angle_w := 360 + angle_w;

  centrex := PaintBox1.Width div 2;
  centrey := PaintBox1.Height div 2;
  for i := 0 to length(obj) - 1 do
    for j := 0 to length(obj[i]) - 1 do
      CalculPoint(obj[i][j], angle_h, angle_v, angle_w, centrex, centrey);

  TestComplet(angle_h, angle_v, angle_w, centrex, centrey);

end;

procedure TForm1.btn_VerifierClick(Sender: TObject);
var
  i, j: integer;
  p, Q: TFace;
  bad: boolean;
begin
  bad := false;
  for i := 0 to length(obj) - 2 do
  begin
    for j := i to length(obj) - 1 do
    begin
      p := obj[i];
      Q := obj[j];
      if ZBufferIntersect(p, Q, PaintBox1.Width div 2, PaintBox1.Height div 2)
        = false then
      begin
        bad := true;
        debug('Bad  = ' + IntToStr(i) + ' / ' + IntToStr(j));
      end;
    end;
  end;
  if bad then
    debug('Model correct');
end;

procedure TForm1.DoChangeAngle(Sender: TObject);
begin
  calculated := false;
  Label1.Caption := IntToStr(scb_h.position);
  Label2.Caption := IntToStr(scb_v.position);
  Label8.Caption := IntToStr(scb_w.position);
  PaintBox1.Invalidate;
end;

procedure TForm1.DoChangeParam(Sender: TObject);
begin
  calculated := false;
end;

procedure TForm1.DoCreate(Sender: TObject);
var
  i: integer;
begin
  DoubleBuffered := true; // pour éviter le flicker
  Label1.Caption := IntToStr(scb_h.position);
  Label2.Caption := IntToStr(scb_v.position);
  Label8.Caption := IntToStr(scb_w.position);
  zoom := StrToInt(edt_Zoom.Text);
  distance := StrToInt(edt_Distance.Text);
  calculated := false;

  // précalcul des sin et cos.
  for i := 0 To 359 do
  begin
    sinus[i] := Sin(i * (pi / 180));
    cosinus[i] := Cos(i * (pi / 180));
  end;
end;

procedure TForm1.DoMemoChange(Sender: TObject);
begin
  if Memo1.Lines.Count > 100 then Memo1.Clear;
end;

procedure TForm1.DoMouseDown(Sender: TObject; Button: TMouseButton;
  Shift: TShiftState; X, Y: Integer);
begin
  MouseIsDown := true;
  PDown := Point(x,y);
  PActually := Point(x,y);
end;

procedure TForm1.DoMouseMove(Sender: TObject; Shift: TShiftState; X,
  Y: Integer);
const
  minimove = 0;
var
  movex, movey : integer;
begin
  if not(MouseIsDown) then  Exit;
  if (abs(PActually.X - X)< minimove) and (abs(PActually.Y - Y)< minimove) then Exit;

  movex := X - PActually.X;
  movey := Y - PActually.Y;
  if (scb_v.Position>-90) and (scb_v.Position<90) then
  scb_h.Position := scb_h.Position - movex
  else
  scb_h.Position := scb_h.Position + movex;

  scb_v.Position := scb_v.Position + movey;
  PActually := Point(x, y);

end;

procedure TForm1.DoMouseUp(Sender: TObject; Button: TMouseButton;
  Shift: TShiftState; x, y: integer);
var
  i, j: integer;
  ptab: array of TPoint;
  r: HRGN;
begin
  // préparation des tableaux de points
  MouseIsDown := false;
  for i := 0 to length(obj) - 1 do
  begin
    SetLength(ptab, length(obj[i]));
    for j := 0 to length(obj[i]) - 1 do
    begin
      ptab[j] := point(obj[i][j].x2d, obj[i][j].y2d);
      r := CreatePolygonRgn(ptab[0], length(obj[i]), WINDING);
    end;
    if PtInRegion(r, x, y) then
    begin
      debug('Face = ' + IntToStr(i));
    end;
    DeleteObject(r);
  end;
end;

procedure TForm1.DoOrthoChange(Sender: TObject);
begin
    PaintBox1.Invalidate;

end;

procedure TForm1.DoWheelDown(Sender: TObject; Shift: TShiftState;
  MousePos: TPoint; var Handled: Boolean);
  //zoom arière sur mouse wheel down
begin
  if ssAlt in Shift then
  edt_Distance.Text := IntToStr(round(StrToInt(edt_Distance.Text)/1.1))
  else
  edt_Zoom.Text := IntToStr(round(StrToInt(edt_Zoom.Text)*1.1));
  btn_AffichageClick(nil);
end;

procedure TForm1.DoWheelUp(Sender: TObject; Shift: TShiftState;
  MousePos: TPoint; var Handled: Boolean);
    //zoom aven sur mouse wheel up
begin
  if ssAlt in Shift then
  edt_Distance.Text := IntToStr(round(StrToInt(edt_Distance.Text)*1.1))
  else
  edt_Zoom.Text := IntToStr(round(StrToInt(edt_Zoom.Text)/1.1));
  btn_AffichageClick(nil);
end;

end.
