
function positivMod(n, mod) {
  return ((n % mod) + mod) % mod;
}

function keepNDec(x, nDec) {

  return Math.round(x * 10 ** nDec) / 10 ** nDec;
}

class point {
  constructor(x, y) {
  
    this.x = x;
    this.y = y;
  }

  copy() {
    return new point(this.x, this.y);
  }

  sameAbsciss(other) {

    return this.x === other.x;
  }

  sameOrdinate(other) {
   
    return this.y === other.y;
  }

  equal(other) {

    return this.sameAbsciss(other) && this.sameOrdinate(other);
  }

  distance(other) {

    return Math.sqrt((this.x - other.y) ** 2 + (this.y - other.y) ** 2);
  }

  addVector(v) {

    let res = new point(this.x + v.x, this.y + v.y);
    return res;
  }

  translate(v) {

    this.x += v.x;
    this.x = keepNDec(this.x, 10);
    this.y += v.y;
    this.y = keepNDec(this.y, 10);
  }
}

class vector {
  constructor(M, N) {
 
    if (M instanceof point) {
      this.x = N.x - M.x;
      this.y = N.y - M.y;
    } else {
      
      this.x = M;
      this.y = N;
    }
  }

  is0() {

    return this.x === 0 && this.y === 0;
  }

  sum(other) {

    let res = new vector(this.x + other.x, this.y + other.y);
    return res;
  }

  product(lambda) {

    let res = new vector(lambda * this.x, lambda * this.y);
    return res;
  }

  scalarProduct(other) {

    return this.x * other.x + this.y * other.y;
  }

  norm() {

    return Math.sqrt(this.scalarProduct(this));
  }

  orthogonalVector() {

    let res = new vector(-this.y, this.x);
    return res;
  }

  polarCoordinate() {

    return [this.norm(), Math.atan(this.y / this.x)];
  }
}

class straightLine {
  constructor(point1, element) {
    if (element instanceof point) {
      if (!point1.equal(element)) {
        this.point1 = point1;
        this.point2 = element;
      } else {
        throw "A straight line can't be definied by two identical points";
      }
    } else {
      if (!element.is0()) {
        this.point1 = point1;
        this.point2 = point1.addVector(element);
      } else {
        throw "A straight line can't be definied by a null vector";
      }
    }
  }

  equation() {

    if (this.point1.sameAbsciss(this.point2)) {
      return [1, 0, -this.point1.x];
    } else {
      let direction =
        (this.point1.y - this.point2.y) / (this.point1.x - this.point2.x);
      let ordinateOrigin = this.point1.y - direction * this.point1.x;
      return [-direction, 1, -ordinateOrigin];
    }
  }

  containPoint(point) {
    let equationLine = this.equation();
    return (
      keepNDec(
        equationLine[0] * point.x + equationLine[1] * point.y + equationLine[2],
        10
      ) === 0
    );
  }
}

class segment {
  constructor(point1, point2) {
    this.point1 = point1;
    this.point2 = point2;
  }

  center() {

    let res = new point(
      (this.point1.x + this.point2.x) / 2,
      (this.point1.y + this.point2.y) / 2
    );
    return res;
  }

  containPoint(point) {
    let vector1 = new vector(this.point1, this.point2);
    let vector2 = new vector(this.point1, point);

    let scalarProd1 = vector1.scalarProduct(vector2);
    let scalarProd2 = vector1.scalarProduct(vector1);

    return scalarProd1 >= 0 && scalarProd1 <= scalarProd2;
  }

  intersect(other) {
    let segmentVector = new vector(this.point1, this.point2);
    let vector1 = new vector(this.point1, other.point1);
    let vector2 = new vector(this.point1, other.point2);
    let vector3 = new vector(this.point2, other.point1);
    let vector4 = new vector(this.point2, other.point2);

    let scalarProd1 = segmentVector.scalarProduct(vector1);
    let scalarProd2 = segmentVector.scalarProduct(vector2);
    let scalarProd3 = segmentVector.scalarProduct(vector3);
    let scalarProd4 = segmentVector.scalarProduct(vector4);

    return !(
      (scalarProd1 < 0 && scalarProd2 < 0) ||
      (scalarProd3 > 0 && scalarProd4 > 0)
    );
  }
}

class polygon {
  constructor(vertices) {
    this.vertices = vertices;
  }

  copy() {
    let newPolygon = new polygon([]);
    this.vertices.forEach((point) => {
      newPolygon.push(point.copy());
    });
    return newPolygon;
  }

  translate(translationVector) {

    this.vertices.forEach((point) => {
      point.translate(translationVector);
    });
  }

  edges() {

    let edges = [];
    let nbVertices = this.vertices.length;
    if (nbVertices > 2) {
      for (let k = 0; k < nbVertices; k++) {
        let edge = new segment(
          this.vertices[k],
          this.vertices[(k + 1) % nbVertices]
        );
        edges.push(edge);
      }
    } else {
      let edge = new segment(this.vertices[0], this.vertices[1]);
      edges.push(edge);
    }
    return edges;
  }

  isoBarycenter() {

    let barycenterAbscissa = 0;
    let barycenterOrdinate = 0;
    this.vertices.forEach((point) => {
      barycenterAbscissa += point.x;
      barycenterOrdinate += point.y;
    });

    let res = new point(
      (1 / this.vertices.length) * barycenterAbscissa,
      (1 / this.vertices.length) * barycenterOrdinate
    );
    return res;
  }

  static separation(other, edge, barycenter) {
    let otherNbVertices = other.vertices.length;
    let segmentLine = new straightLine(edge.point1, edge.point2);
    let equation = segmentLine.equation();

    let thisSide =
      equation[0] * barycenter.x + equation[1] * barycenter.y + equation[2];

    let pointSideSet = [];
    let pointSide = [];

    let pointOnSepartor = [];

    for (let k = 0; k < otherNbVertices; k++) {
      pointSide =
        equation[0] * other.vertices[k].x +
        equation[1] * other.vertices[k].y +
        equation[2];
      pointSideSet.push(pointSide);
      if (keepNDec(pointSide, 10) === 0) {
        pointOnSepartor.push(other.vertices[k]);
      }
    }

    let commonPoint = false;
    if (pointOnSepartor.length == 1) {
      if (edge.containPoint(pointOnSepartor[0])) {
        commonPoint = true;
      }
    } else if (pointOnSepartor.length == 2) {
      let alignSegment = new segment(pointOnSepartor[0], pointOnSepartor[1]);
      if (edge.intersect(alignSegment)) {
        commonPoint = true;
      }
    }

    if (commonPoint) {
      return false;
    } else {
      let minPointSide = Math.min.apply(Math, pointSideSet);
      let maxPointSide = Math.max.apply(Math, pointSideSet);

      if (keepNDec(thisSide, 10) == 0) {
        return keepNDec(minPointSide, 10) * keepNDec(maxPointSide, 10) >= 0;
      } else {
        return (
          keepNDec(thisSide, 10) * keepNDec(maxPointSide, 10) <= 0 &&
          keepNDec(minPointSide, 10) * keepNDec(thisSide, 10) <= 0
        );
      }
    }
  }

  sat(other) {

    let thisEdges = this.edges();
    let otherEdges = other.edges();

    let thisBarycenter = this.isoBarycenter();
    let otherBarycenter = other.isoBarycenter();

    let isSeparated = false;
    let cpt = 0;
    do {

      isSeparated = polygon.separation(other, thisEdges[cpt], thisBarycenter);
      cpt++;
    } while ((cpt < thisEdges.length) & !isSeparated);

    if (!isSeparated) {

      cpt = 0;
      do {
        isSeparated = polygon.separation(
          this,
          otherEdges[cpt],
          otherBarycenter
        );
        cpt++;
      } while ((cpt < otherEdges.length) & !isSeparated);
    }

    return isSeparated;
  }
}


let gameScore = 0;
let lastScoreTime = 0;

function updateScore() {
  if (heroInstance.hasStarted && !heroInstance.isDead && !heroInstance.haveFinished) {
    if (Date.now() - lastScoreTime >= 1000) { // 1000 ms = 1 segundo
      gameScore += 10;
      lastScoreTime = Date.now();
    }
  }
}

function drawScore() {
  drawingInstance.ctx.font = "60px Orbitron";
  drawingInstance.ctx.shadowColor = "rgba(224, 231, 34)";
  drawingInstance.ctx.shadowBlur = 20;
  drawingInstance.ctx.fillStyle = "rgba(224, 231, 34, 1)";
  drawingInstance.ctx.textAlign = "right";
  drawingInstance.ctx.fillText("Puntos: " + gameScore, drawingInstance.width - 20, 85);
  drawingInstance.ctx.shadowBlur = 0;
}

class square extends polygon {

  constructor(element1, element2) {

    let point1;
    let point2;
    let point3;
    let point4;

    let direction;
    let polarDirection;
    let diagonal;
    let center;

    if (element2 instanceof point) {
      point1 = element1;
      point2 = element2;
      direction = new vector(point1, point2);
      polarDirection = direction.polarCoordinate();
      point3 = point2.addVector(direction.orthogonalVector());
      point4 = point3.addVector(
        direction.orthogonalVector().orthogonalVector()
      );
      diagonal = new segment(point1, point3);
      center = diagonal.center();
    } else {
      polarDirection = element2;
      direction = new vector(
        polarDirection[0] * Math.cos(polarDirection[1]),
        polarDirection[0] * Math.sin(polarDirection[1])
      );

      point1 = new point(0, 0);
      point2 = point1.addVector(direction);
      point3 = point2.addVector(direction.orthogonalVector());
      point4 = point3.addVector(
        direction.orthogonalVector().orthogonalVector()
      );


      diagonal = new segment(point1, point3);
      let initialCenter = diagonal.center();
      center = element1;
      let translationVector = new vector(initialCenter, center);

      point1.translate(translationVector);
      point2.translate(translationVector);
      point3.translate(translationVector);
      point4.translate(translationVector);
    }
    super([point1, point2, point3, point4]);
    this.center = center;
    this.polarDirection = polarDirection;
  }

  copy() {
    let newSquare = new square([
      this.vertices[0].copy(),
      this.vertices[0].copy(),
    ]);
    newSquare.center = this.center.copy();
    newSquare.polarDirection = this.polarDirection.slice();
    return newSquare;
  }

  rotate(angle) {

    this.polarDirection[1] += angle;
    let direction = new vector(
      this.polarDirection[0] * Math.cos(this.polarDirection[1]),
      this.polarDirection[0] * Math.sin(this.polarDirection[1])
    );

    this.vertices[0] = new point(0, 0);
    this.vertices[1] = this.vertices[0].addVector(direction);
    this.vertices[2] = this.vertices[1].addVector(direction.orthogonalVector());
    this.vertices[3] = this.vertices[2].addVector(
      direction.orthogonalVector().orthogonalVector()
    );

    let diagonal = new segment(this.vertices[0], this.vertices[2]);
    let initialCenter = diagonal.center();
    let translationVector = new vector(initialCenter, this.center);

    this.vertices[0].translate(translationVector);
    this.vertices[1].translate(translationVector);
    this.vertices[2].translate(translationVector);
    this.vertices[3].translate(translationVector);
  }

  translate(transactionVector) {

    super.translate(transactionVector);
    this.center.translate(transactionVector);
  }

  getLowestPointIndex() {

    let lowestPoint = new point(Infinity, Infinity);

    for (let k = 0; k < this.vertices.length; k++) {
      if (keepNDec(this.vertices[k].y, 6) < keepNDec(lowestPoint.y, 6)) {

        lowestPoint = this.vertices[k];
      }
    }

    let res = [];
    for (let k = 0; k < this.vertices.length; k++) {
      if (keepNDec(lowestPoint.y, 6) === keepNDec(this.vertices[k].y, 6)) {

        res.push(k);
      }
    }

    if (res.length == 2) {

      if (this.vertices[res[0]].x > this.vertices[res[1]].x) {
        res = [res[1], res[0]];
      }
    }

    return res;
  }
}



class hero {
  constructor(
    positionCenter,
    positionPolarCoordinates,
    vx,
    vy0,
    xJump,
    yJump,
    g,
    t,
    isJumping
  ) {



    let intialPosition = new point(positionCenter[0], positionCenter[1]);
    this.body = new square(intialPosition, positionPolarCoordinates.slice());


    let footPoint = this.body.getLowestPointIndex();
    if (footPoint.length == 2) {

      let footPoint1 = this.body.vertices[footPoint[0]].copy();
      let footPoint2 = this.body.vertices[footPoint[1]].copy();
      let foot1 = new polygon([footPoint1, footPoint2]);
      this.foot = [foot1];
    } else {

      let footPoint1 = this.body.vertices[footPoint[0]].copy();
      let footPoint2 =
        this.body.vertices[positivMod(footPoint[0] + 1, 4)].copy();
      let footPoint3 =
        this.body.vertices[positivMod(footPoint[0] - 1, 4)].copy();
      let foot1 = new polygon([footPoint1, footPoint2]);
      let foot2 = new polygon([footPoint1.copy(), footPoint3]);
      this.foot = [foot1, foot2];
    }



    this.vx = vx;
    this.vy0 = vy0;
    this.xJump = xJump;
    this.yJump = yJump;
    this.g = g;
    this.t = t;
    this.isJumping = isJumping;

    this.startJumpPosition = this.body.center.copy();



    this.hasStarted = false;
    this.isDead = false;
    this.haveFinished = false;



    this.deathParticle = [];
  }

  rotate(angle) {

    this.body.rotate(angle);

    let footPoint = this.body.getLowestPointIndex();
    if (footPoint.length == 2) {
      let footPoint1 = this.body.vertices[footPoint[0]].copy();
      let footPoint2 = this.body.vertices[footPoint[1]].copy();
      let foot1 = new polygon([footPoint1, footPoint2]);
      this.foot = [foot1];
    } else {
      let footPoint1 = this.body.vertices[footPoint[0]].copy();
      let footPoint2 =
        this.body.vertices[positivMod(footPoint[0] + 1, 4)].copy();
      let footPoint3 =
        this.body.vertices[positivMod(footPoint[0] - 1, 4)].copy();
      let foot1 = new polygon([footPoint1, footPoint2]);
      let foot2 = new polygon([footPoint1.copy(), footPoint3]);
      this.foot = [foot1, foot2];
    }
  }

  translate(transactionVector) {

    this.body.translate(transactionVector);
    this.foot.forEach((footValue) => {
      footValue.translate(transactionVector);
    });
  }

  footContactWithRoof(previousFoot, platformInstance) {


    let cpt = 0;

    let footPolygon;
    for (let k = 0; k < this.foot.length; k++) {
      let lineTest = new straightLine(
        this.foot[k].vertices[1],
        previousFoot[k][0]
      );

      if (lineTest.containPoint(previousFoot[k][1])) {
        footPolygon = new polygon([
          this.foot[k].vertices[1],
          previousFoot[k][0],
        ]);

      } else {
        footPolygon = new polygon([
          this.foot[k].vertices[0],
          this.foot[k].vertices[1],
          previousFoot[k][1],
          previousFoot[k][0],
        ]);
      }

      if (!footPolygon.sat(platformInstance.roof)) {

        cpt++;
      }
    }

    return cpt > 0;
  }

  move(gridInstance) {


    let previousFoot = [];
    this.foot.forEach((foot) => {
      previousFoot.push([foot.vertices[0].copy(), foot.vertices[1].copy()]);
    });

    let tSauv = this.t;
    let yPosSauv = this.body.center.y;
    let xPosSauv = this.body.center.x;


    let dt = frameTimeDiff.dt;


    let translationVector = new vector(
      this.vx * dt,
      Math.max((-1 / 2) * this.g * dt * (2 * this.t + dt) + this.vy0 * dt, -1)
    );


    this.translate(translationVector);


    this.t += dt;

    let deadContactElementCenter = [];

    let floorContactElementCenter = [];


    let aroundGrid = gridInstance.grid.slice(
      Math.max(Math.floor(this.body.center.x - 1), 0),
      Math.floor(this.body.center.x + 2)
    );

    aroundGrid.forEach((col) => {

      if (col != undefined) {

        col.forEach((element) => {
          if (element instanceof platform) {
            if (!this.body.sat(element.platform)) {

              if (this.footContactWithRoof(previousFoot, element)) {
                floorContactElementCenter.push(element.platform.center.y);
              } else {
                deadContactElementCenter.push(element.platform.center.y);
              }
            }
          } else if (element instanceof peak) {

            if (!this.body.sat(element.peak)) {
              deadContactElementCenter.push(element.center.y);
            }
          } else if (element instanceof ending) {

            if (!this.body.sat(element.ending)) {
              this.haveFinished = true;
            }
          }
        });
      }
    });

    let maxDeadContactCenter, maxFloorContactCenter;

    if (deadContactElementCenter.length > 0) {
      maxDeadContactCenter = deadContactElementCenter.reduce(function (a, b) {
        return Math.max(a, b);
      });
    } else {
      maxDeadContactCenter = -Infinity;
    }

    if (floorContactElementCenter.length > 0) {
      maxFloorContactCenter = floorContactElementCenter.reduce(function (a, b) {
        return Math.max(a, b);
      });
    } else {
      maxFloorContactCenter = -Infinity;
    }

    if (deadContactElementCenter.length > 0) {
      if (floorContactElementCenter.length === 0) {
        this.isDead = true;
      } else {
        if (maxDeadContactCenter > maxFloorContactCenter) {
          this.isDead = true;
        }
      }
    }

    if (this.body.center.y < 0) {

      this.isDead = true;
      let newCenter = new point(this.body.center.x, 0.5);
      let translateVector = new vector(this.body.center, newCenter);
      this.translate(translateVector);
    }

    if (!this.isDead) {
      if (floorContactElementCenter.length > 0) {

        let newXPosition;
        if (this.isJumping) {

          let a = -this.g / 2;
          let b = this.vy0 - this.g * tSauv;
          let c = -(maxFloorContactCenter + 1 - yPosSauv);
          let delta = b ** 2 - 4 * a * c;
          let newDt = Math.max(
            (-b - Math.sqrt(delta)) / (2 * a),
            (-b + Math.sqrt(delta)) / (2 * a)
          );

          newXPosition = xPosSauv + newDt * this.vx;
        } else {

          newXPosition = this.body.center.x;
        }

        let newCenter = new point(newXPosition, maxFloorContactCenter + 1);
        let translateVector = new vector(this.body.center, newCenter);
        this.translate(translateVector);


        this.rotate(2 * Math.PI - this.body.polarDirection[1]);


        this.g = 0; // g is compensated by newton 3rd law
        this.vy0 = 0; // the vertical movment stop
        this.t = 0; // t is not mandatory in equation anymore, and will be 0 at the begning of a jump or a fall
        this.isJumping = false; // the gamer can jump
      } else {
        this.g = (2 * this.yJump) / (this.xJump / (2 * this.vx)) ** 2; // no compensation by newton 3rd law
        this.rotate(
          -Math.PI / (((1 / frameTimeDiff.dt) * this.xJump) / (2 * this.vx) + 2)
        ); // to look pretty : the hero rotate when not on a roof
        this.isJumping = true; // can't jump before the end of the jump /fall
      }
    }
  }

  jump() {
    /* Modify physical constant in order to the next move (methode move) is a jump */
    if (!this.isJumping) {
      this.isJumping = true; // during a jump, the hero can't jump anymore
      // can't jump if already jumping
      this.g = (2 * this.yJump) / (this.xJump / (2 * this.vx)) ** 2; // no compensation by newton 3rd law
      this.vy0 = (2 * this.yJump) / (this.xJump / (2 * this.vx)); // to look pretty : the hero rotate when not on a roof
      this.t = 0; // at the begning of a jump, t must be = 0 in order to the equation are ok
      this.startJumpPosition = this.body.center.copy();
    }
  }

  setDeathParticle() {
    /*  When death, create explosion particles */
    for (let k = 0; k < 40; k++) {
      this.deathParticle.push({
        position: this.body.center.copy(),
        angle: 2 * Math.PI * Math.random(),
        maxProjection: 2 * Math.random(),
      });
    }
  }
}


class platform {

  constructor(x, y) {
    /* x,y are the coordinates of the lowest leftest point of the platform */
    this.col = Math.floor(x);
    let platformCenter = new point(x + 1 / 2, y + 1 / 2);
    this.platform = new square(platformCenter, [1, 0]);
    this.roof = new polygon([
      this.platform.vertices[2],
      this.platform.vertices[3],
    ]);
  }
}

class peak {

  constructor(x, y, orientation) {
    /* x, y are the coordinates of the lowest leftest point of the square in which the triangle is inscribed */
    let point1, point2, point3;
    switch (orientation) {
      case "up":
        point1 = new point(x, y);
        point2 = new point(x + 1, y);
        point3 = new point(x + 1 / 2, y + 1);
        break;
      case "down":
        point1 = new point(x, y + 1);
        point2 = new point(x + 1, y + 1);
        point3 = new point(x + 1 / 2, y);
        break;
      case "left":
        point1 = new point(x + 1, y + 1);
        point2 = new point(x + 1, y);
        point3 = new point(x, y + 1 / 2);
        break;
      case "right":
        point1 = new point(x, y + 1);
        point2 = new point(x, y);
        point3 = new point(x + 1, y + 1 / 2);
        break;
    }
    this.col = Math.floor(x);
    this.peak = new polygon([point1, point2, point3]);
    this.center = new point(x + 1 / 2, y + 1 / 2);
  }
}

class ending {

  constructor(position) {

    let point1 = new point(position, 0);
    let point2 = new point(position + 1, 0);
    let point3 = new point(position + 1, 10);
    let point4 = new point(position, 10);

    this.ending = new polygon([point1, point2, point3, point4]);
    this.col = Math.floor(position);
  }
}

class checkPoint {

  constructor(heroInstance) {
    /*  the position of the checkpoint depends of the position of the hero */
    this.x = heroInstance.body.center.x;
    this.y = heroInstance.body.center.y;
    this.col = Math.floor(this.x);
  }

  update(heroInstance) {
    this.x = heroInstance.body.center.x;
    this.y = heroInstance.body.center.y;
    this.col = Math.floor(this.x);
  }
}

class grid {

  constructor() {
    /* by default a grid is empty, we add element using methodes */
    this.grid = [];
  }

  addPlatform(x, y) {
    let platformInstance = new platform(x, y);
    /*  Given a platform, place a platform in the grid */
    if (this.grid[platformInstance.col] != undefined) {
      this.grid[platformInstance.col].push(platformInstance);
    } else {
      this.grid[platformInstance.col] = [platformInstance];
    }

    return [x, y];
  }

  addPeak(x, y, direction) {
    let peakInstance = new peak(x, y, direction);
    /*  Given a peak, place a platform in the grid */
    if (this.grid[peakInstance.col] != undefined) {
      this.grid[peakInstance.col].push(peakInstance);
    } else {
      this.grid[peakInstance.col] = [peakInstance];
    }

    return [x, y];
  }

  addEnding(endingInstance) {
    /*  Given a ending, place a platform in the grid */
    if (this.grid[endingInstance.col] != undefined) {
      this.grid[endingInstance.col].push(endingInstance);
    } else {
      this.grid[endingInstance.col] = [endingInstance];
    }
  }

  addCheckPoint(checkPointInstance) {
    /*  Given a checkPoint, place a platform in the grid */
    if (this.grid[checkPointInstance.col] != undefined) {
      this.grid[checkPointInstance.col].push(checkPointInstance);
    } else {
      this.grid[checkPointInstance.col] = [checkPointInstance];
    }
  }

  removeCol(start, end) {
    /* remove the entire element of grid from col start to col end exclude (begning by 0) */
    for (let k = start; k < end; k++) {
      this.grid[k] = undefined;
    }
  }

  defaultGrid(size) {
    /* create a default frid with a default ground of platform on first row. */
    this.grid = [];
    for (let k = 4; k < size; k++) {
      this.addPlatform(k, 4);
    }
  }
}



class backGround {

  constructor() {
    this.city = document.getElementById("city");
    this.cityReverse = document.getElementById("city-reverse");

    this.t = 0;

  }
}

class drawing {


  constructor(heroInstance) {
    // principal canvas
    let canvas = document.getElementById("canvas-game");
    this.ctx = canvas.getContext("2d");
    this.width = document.getElementById("game-interface").offsetWidth;
    this.height = document.getElementById("game-interface").offsetHeight;
    this.ctx.canvas.width = this.width;
    this.ctx.canvas.height = this.height;
    this.unity = this.width / 40; // unity of the game. for exemple, segment (0,0), (1,0) will have unity length when draw
    this.heroCenterXPosition = 0; // use to make the grid scroll with the hero on x. update in method setGridPosition
    this.heroAjustYPosition = 0; // use to make the grid scrill with the hero on y. update in method setGridPosition
    this.deathAnimationTime = 0.3;
    this.winAnimationTime = 2;

    // background canvas
    let canvasBack = document.getElementById("canvas-background");
    this.ctxBack = canvasBack.getContext("2d");
    this.ctxBack.canvas.width = this.width;
    this.ctxBack.canvas.height = this.height;
    this.backGroundTimeScroll = 0;
    this.backgroundSpeed = (heroInstance.vx / 6) * this.unity;

    this.backgroundImageCity = document.getElementById("city");
    this.backgroundImageCityReverse = document.getElementById("city-reverse");

    //
    this.sound = {
      backGroundMusic: document.getElementById("backGroundMusic"),
      checkpointSaveTime: 0,
      deathSound: document.getElementById("deathSound"),
    };
  }

  setGridPosition(heroInstance) {

    this.heroCenterXPosition = heroInstance.body.center.x - 10;
    if (heroInstance.body.center.y + this.heroAjustYPosition < 5) {
      this.heroAjustYPosition = Math.min(5 - heroInstance.body.center.y, 0);
    } else if (
      this.height / this.unity -
      heroInstance.body.center.y -
      this.heroAjustYPosition <
      5
    ) {
      this.heroAjustYPosition =
        this.height / this.unity - heroInstance.body.center.y - 5;
    }
  }

  gridAbscissa(x) {
    
    return x * this.unity;
  }

  gridOrdinate(y) {

    return this.height - y * this.unity;
  }

  drawSquareNeonStyle(elementToDraw, r, g, b) {
    // Primer color (borde)
    this.ctx.strokeStyle = `rgba(0,0,0,1)`;
    this.ctx.lineWidth = 4;
    this.ctx.stroke(elementToDraw);

    // Segundo color (relleno, diferente)
    // Genera un color complementario simple para el relleno
    const r2 = 255 - r;
    const g2 = 255 - g;
    const b2 = 255 - b;
    this.ctx.fillStyle = `rgba(${r2},${g2},${b2},0.6)`;
    this.ctx.fill(elementToDraw);
  }

  drawTextNeonStyle(
    text,
    r,
    g,
    b,
    alpha,
    fontSize,
    position,
    align = "center"
  ) {
    this.ctx.shadowColor = "rgba(" + r + "," + g + "," + b + ")";
    this.ctx.shadowBlur = 10;
    this.ctx.font = fontSize + "px Orbitron";
    this.ctx.fillStyle = "rgba(" + r + "," + g + "," + b + "," + alpha + ")";
    this.ctx.textAlign = align;
    this.ctx.fillText(
      text,
      this.gridAbscissa(position[0] - this.heroCenterXPosition),
      this.gridOrdinate(position[1] + this.heroAjustYPosition)
    );
    this.ctx.shadowBlur = 0;
  }

  drawMovingtext(
    text,
    r,
    g,
    b,
    alpha,
    fontSize,
    position1,
    position2,
    heroInstance
  ) {
    if (
      (position1[0] <= heroInstance.body.center.x + 40) &
      (position1[0] > heroInstance.body.center.x - 5)
    ) {
      this.drawTextNeonStyle(text, r, g, b, alpha, fontSize, position1, "left");
    } else if (
      position2[0] > heroInstance.body.center.x &&
      position1[0] < heroInstance.body.center.x
    ) {
      this.drawTextNeonStyle(
        text,
        r,
        g,
        b,
        alpha,
        fontSize,
        [heroInstance.body.center.x - 5, position1[1]],
        "left"
      );
    } else if (position2[0] > heroInstance.body.center.x - 40) {
      this.drawTextNeonStyle(
        text,
        r,
        g,
        b,
        alpha,
        fontSize,
        [position2[0] - 5, position2[1]],
        "left"
      );
    }
  }

  drawHero(heroInstance) {

    let heroBody = new Path2D();
    heroBody.moveTo(
      this.gridAbscissa(
        heroInstance.body.vertices[0].x - this.heroCenterXPosition
      ),
      this.gridOrdinate(
        heroInstance.body.vertices[0].y + this.heroAjustYPosition
      )
    );
    heroBody.lineTo(
      this.gridAbscissa(
        heroInstance.body.vertices[1].x - this.heroCenterXPosition
      ),
      this.gridOrdinate(
        heroInstance.body.vertices[1].y + this.heroAjustYPosition
      )
    );
    heroBody.lineTo(
      this.gridAbscissa(
        heroInstance.body.vertices[2].x - this.heroCenterXPosition
      ),
      this.gridOrdinate(
        heroInstance.body.vertices[2].y + this.heroAjustYPosition
      )
    );
    heroBody.lineTo(
      this.gridAbscissa(
        heroInstance.body.vertices[3].x - this.heroCenterXPosition
      ),
      this.gridOrdinate(
        heroInstance.body.vertices[3].y + this.heroAjustYPosition
      )
    );
    heroBody.closePath();
    this.drawSquareNeonStyle(heroBody, 254, 1, 154);
  }

  drawGrid(gridInstance, heroInstance) {

    let minColToKeep = Math.floor(Math.max(this.heroCenterXPosition, 0));
    let maxColToKeep = Math.floor(heroInstance.body.center.x + 40);

    let platformDraw = new Path2D();
    let peakDraw = new Path2D();
    let checkPointDraw = new Path2D();
    let endingDraw = new Path2D();
    let endingPosition;

    gridInstance.grid.slice(minColToKeep, maxColToKeep).forEach((gridRow) => {

      if (gridRow != undefined) {
        gridRow.forEach((element) => {
          if (element instanceof platform) {
            platformDraw.rect(
              this.gridAbscissa(
                element.platform.vertices[3].x - this.heroCenterXPosition
              ),
              this.gridOrdinate(
                element.platform.vertices[3].y + this.heroAjustYPosition
              ),
              this.unity,
              this.unity
            );

            platformDraw.closePath();
          } else if (element instanceof peak) {
            peakDraw.moveTo(
              this.gridAbscissa(
                element.peak.vertices[0].x - this.heroCenterXPosition
              ),
              this.gridOrdinate(
                element.peak.vertices[0].y + this.heroAjustYPosition
              )
            );
            peakDraw.lineTo(
              this.gridAbscissa(
                element.peak.vertices[1].x - this.heroCenterXPosition
              ),
              this.gridOrdinate(
                element.peak.vertices[1].y + this.heroAjustYPosition
              )
            );
            peakDraw.lineTo(
              this.gridAbscissa(
                element.peak.vertices[2].x - this.heroCenterXPosition
              ),
              this.gridOrdinate(
                element.peak.vertices[2].y + this.heroAjustYPosition
              )
            );
            peakDraw.closePath();
          } else if (element instanceof checkPoint) {
            checkPointDraw.moveTo(
              this.gridAbscissa(element.x - this.heroCenterXPosition),
              this.gridOrdinate(element.y + this.heroAjustYPosition)
            );
            checkPointDraw.lineTo(
              this.gridAbscissa(element.x - this.heroCenterXPosition),
              this.gridOrdinate(element.y + 2 + this.heroAjustYPosition)
            );
            checkPointDraw.lineTo(
              this.gridAbscissa(element.x + 1 - this.heroCenterXPosition),
              this.gridOrdinate(element.y + 3 / 2 + this.heroAjustYPosition)
            );
            checkPointDraw.lineTo(
              this.gridAbscissa(element.x - this.heroCenterXPosition),
              this.gridOrdinate(element.y + 1 + this.heroAjustYPosition)
            );
            checkPointDraw.closePath();
          } else if (element instanceof ending) {
            endingDraw.moveTo(
              this.gridAbscissa(element.col - 1 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(14 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col + 3 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(14 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col + 1 - this.heroCenterXPosition),
              this.gridOrdinate(6 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col + 3 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(6 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col + 1 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(5 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col - 1 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(6 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col - this.heroCenterXPosition),
              this.gridOrdinate(6 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col - 1 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(14 + this.heroAjustYPosition)
            );
            endingDraw.lineTo(
              this.gridAbscissa(element.col + 3 / 2 - this.heroCenterXPosition),
              this.gridOrdinate(14 + this.heroAjustYPosition)
            );

            checkPointDraw.closePath();
            endingPosition = element.col;
          }
        });
      }
    });

    // print elements
    this.drawSquareNeonStyle(peakDraw, 255, 7, 58);
    this.drawSquareNeonStyle(platformDraw, 255, 254, 242);
    this.drawSquareNeonStyle(checkPointDraw, 224, 231, 34);
    this.drawMovingtext(
      "Presione ESPACIO para comenzar",
      70,
      102,
      255,
      1,
      60,
      [0, 16],
      [10, 16],
      heroInstance
    );
    this.drawMovingtext(
      "Presione ESPACIO para saltar",
      70,
      102,
      255,
      1,
      60,
      [50, 16],
      [80, 16],
      heroInstance
    );
    this.drawMovingtext(
      "Mantenga ESPACIO para saltar largo  ",
      70,
      102,
      255,
      1,
      60,
      [110, 16],
      [140, 16],
      heroInstance
    );
    this.drawMovingtext(
      "Presione S para un punto de control",
      70,
      102,
      255,
      1,
      60,
      [170, 16],
      [200, 16],
      heroInstance
    );

    if (
      Date.now() - frameTimeDiff.endingBegin <
      drawingInstance.winAnimationTime * 1500 ||
      (!heroInstance.isDead && !heroInstance.haveFinished)
    ) {

      this.drawSquareNeonStyle(endingDraw, 57, 255, 20);
      if (endingPosition != undefined) {
        this.drawTextNeonStyle("GANASTE!", 57, 255, 20, 1, 80, [
          endingPosition + 1 / 2,
          16,
        ]);
        if (checkPointCounter > 0) {
          this.drawTextNeonStyle(
            "¿PUEDES HACERLO SIN PUNTOS DE CONTROL?",
            224,
            231,
            34,
            1,
            20,
            [endingPosition + 4 / 2, 15.5]
          );
        }
      }
    }
  }

  drawCheckpointCounter() {

    this.ctxBack.font = "80px Orbitron";
    this.ctxBack.shadowColor = "rgba(224,231,34)";
    this.ctxBack.shadowBlur = 20;
    this.ctxBack.fillStyle = "rgba(224,231,34,1)";
    this.ctxBack.fillText(checkPointCounter, 15, 85);
    this.ctxBack.shadowBlur = 0;
  }

  drawPressToRestart(t) {
    this.ctx.font = "40px Orbitron";
    this.ctx.fillStyle = "rgba(57, 255, 20," + (Math.sin(t) / 3 + 2 / 3) + ")";
    this.ctx.textAlign = "center";
    this.ctx.fillText(
      "Presione espacio para comenzar de nuevo",
      this.ctx.canvas.width / 2,
      this.ctx.canvas.height / 2
    );
  }

  deathAnimation(heroInstance) {

    let translationVector;
    let explosion = new Path2D();

    heroInstance.deathParticle.forEach((particle) => {
      translationVector = new vector(
        Math.cos(particle.angle),
        Math.sin(particle.angle)
      );
      translationVector = translationVector.product(
        particle.maxProjection / (this.deathAnimationTime / frameTimeDiff.dt)
      );
      particle.position.translate(translationVector);
      explosion.moveTo(
        this.gridAbscissa(particle.position.x - this.heroCenterXPosition) + 2,
        this.gridOrdinate(particle.position.y + this.heroAjustYPosition)
      );
      explosion.arc(
        this.gridAbscissa(particle.position.x - this.heroCenterXPosition),
        this.gridOrdinate(particle.position.y + this.heroAjustYPosition),
        2,
        0,
        2 * Math.PI
      );
      explosion.closePath();
    });

    this.drawSquareNeonStyle(explosion, 254, 1, 154);
  }

  drawBackGround() {
    let ratio =
      this.backgroundImageCity.width / this.backgroundImageCity.height;
    let imageHeight = this.ctxBack.canvas.height;
    let imageWidth = imageHeight * ratio;

    this.backGroundTimeScroll += frameTimeDiff.dt * this.backgroundSpeed;

    let xPosition = Math.floor(this.backGroundTimeScroll / imageWidth);

    this.ctxBack.clearRect(
      0,
      0,
      this.ctxBack.canvas.width,
      this.ctxBack.canvas.height
    );

    if (xPosition % 2 == 0) {

      this.ctxBack.drawImage(
        this.backgroundImageCity,
        xPosition * imageWidth - this.backGroundTimeScroll,
        0,
        imageWidth,
        imageHeight
      );
      this.ctxBack.drawImage(
        this.backgroundImageCityReverse,
        (xPosition + 1) * imageWidth - this.backGroundTimeScroll,
        0,
        imageWidth,
        imageHeight
      );
      this.ctxBack.drawImage(
        this.backgroundImageCity,
        (xPosition + 2) * imageWidth - this.backGroundTimeScroll,
        0,
        imageWidth,
        imageHeight
      );
    } else {
      this.ctxBack.drawImage(
        this.backgroundImageCityReverse,
        xPosition * imageWidth - this.backGroundTimeScroll,
        0,
        imageWidth,
        imageHeight
      );
      this.ctxBack.drawImage(
        this.backgroundImageCity,
        (xPosition + 1) * imageWidth - this.backGroundTimeScroll,
        0,
        imageWidth,
        imageHeight
      );
      this.ctxBack.drawImage(
        this.backgroundImageCityReverse,
        (xPosition + 2) * imageWidth - this.backGroundTimeScroll,
        0,
        imageWidth,
        imageHeight
      );
    }
  }
}



function restart(parameters) {

  heroInstance = new hero(
    parameters[0],
    parameters[1],
    parameters[2],
    parameters[3],
    parameters[4],
    parameters[5],
    parameters[6],
    parameters[7],
    parameters[8]
  );
  drawingInstance.backGroundTimeScroll = backGroundPositionSauv.save;
  drawingInstance.ctx.clearRect(
    0,
    0,
    drawingInstance.width,
    drawingInstance.height
  );
  drawingInstance.drawBackGround();
  drawingInstance.setGridPosition(heroInstance);
  drawingInstance.drawHero(heroInstance);
  drawingInstance.drawGrid(gridInstance, heroInstance);
}

function deathFinish() {

  drawingInstance.sound.deathSound.play();
  drawingInstance.ctx.clearRect(
    0,
    0,
    drawingInstance.width,
    drawingInstance.height
  );
  drawingInstance.setGridPosition(heroInstance);
  drawingInstance.drawGrid(gridInstance, heroInstance);
  drawingInstance.deathAnimation(heroInstance);
  if (
    Date.now() - frameTimeDiff.endingBegin <
    drawingInstance.deathAnimationTime * 1000
  ) {
    requestAnimationFrame(deathFinish);
  } else {
    if (checkPointCounter > 0) {

      restart(gameParameters.saved);
      heroInstance.isDead = false;
      heroInstance.hasStarted = true;
      heroInstance.havefinished = false;
      drawingInstance.sound.backGroundMusic.pause();
      drawingInstance.sound.backGroundMusic.currentTime =
        drawingInstance.sound.checkpointSaveTime;
      drawingInstance.sound.backGroundMusic.play();
      frameTimeDiff.lastTime = Date.now();
      game();
    } else {
      drawingInstance.sound.backGroundMusic.pause();
      drawingInstance.sound.backGroundMusic.currentTime = 0;
      frameTimeDiff.startBegin = Date.now();
      backGroundPositionSauv.save = 0;
      level1(gridInstance, heroInstance);
      restart(gameParameters.initial);
    }
  }
}

function winFinish() {

  heroInstance.move(gridInstance);
  drawingInstance.ctx.clearRect(
    0,
    0,
    drawingInstance.width,
    drawingInstance.height
  );
  drawingInstance.setGridPosition(heroInstance);
  drawingInstance.drawHero(heroInstance);
  drawingInstance.drawGrid(gridInstance, heroInstance);
  if (
    Date.now() - frameTimeDiff.endingBegin <
    drawingInstance.winAnimationTime * 1500
  ) {

    frameTimeDiff.dt *= 0.95;
    requestAnimationFrame(winFinish);
  } else {

    if (keys.Space) {
      drawingInstance.sound.backGroundMusic.pause();
      drawingInstance.sound.backGroundMusic.currentTime = 0;
      frameTimeDiff.startBegin = Date.now();
      backGroundPositionSauv.save = 0;
      level1(gridInstance, heroInstance);
      restart(gameParameters.initial);
    } else {
      requestAnimationFrame(winFinish);
      drawingInstance.drawPressToRestart(
        (Math.PI * 2 * (Date.now() - frameTimeDiff.endingBegin)) / 1000
      );
    }
  }
}

function game() {

  frameTimeDiff.dt = (Date.now() - frameTimeDiff.lastTime) / 1000;
  frameTimeDiff.lastTime = Date.now();


  if (!heroInstance.isDead && !heroInstance.haveFinished) {
    if (keys.Space && Date.now() - frameTimeDiff.startBegin > 500) {
      heroInstance.jump();
    }
    if (keys.KeyS && Date.now() - frameTimeDiff.lastCheckPoint > 200) {
      addCheckPoint();
      checkPointCounter++;
      checkPointValue.checkpoint.update(heroInstance);
      gridInstance.addCheckPoint(checkPointValue.checkpoint);
      frameTimeDiff.lastCheckPoint = Date.now();
      backGroundPositionSauv.save = drawingInstance.backGroundTimeScroll;
      drawingInstance.sound.checkpointSaveTime =
        drawingInstance.sound.backGroundMusic.currentTime;
    }
    heroInstance.move(gridInstance);
    drawingInstance.ctx.clearRect(
      0,
      0,
      drawingInstance.width,
      drawingInstance.height
    );
    drawingInstance.drawBackGround();
    drawingInstance.setGridPosition(heroInstance);
    drawingInstance.drawHero(heroInstance);
    drawingInstance.drawGrid(gridInstance, heroInstance);
    drawingInstance.drawCheckpointCounter();
    requestAnimationFrame(game);
  } else {
    frameTimeDiff.endingBegin = Date.now();
    if (heroInstance.isDead) {
      heroInstance.setDeathParticle();
      deathFinish();
    } else {
      winFinish();
    }
  }
}


function addCheckPoint() {

  gameParameters.saved = [
    [heroInstance.body.center.x, heroInstance.body.center.y],
    heroInstance.body.polarDirection.slice(),
    heroInstance.vx,
    heroInstance.vy0,
    heroInstance.xJump,
    heroInstance.yJump,
    heroInstance.g,
    heroInstance.t,
    heroInstance.isJumping,
  ];
}

function keyEventHandler(event) {

  keys[event.code] = event.type === "keydown";
  event.preventDefault();
  if (!heroInstance.hasStarted && keys.Space) {
    drawingInstance.sound.backGroundMusic.play();
    heroInstance.hasStarted = true;
    checkPointCounter = 0;
    gameParameters.saved = [
      gameParameters.initial[0].slice(),
      gameParameters.initial[1].slice(),
      gameParameters.initial.slice(2),
    ];
    frameTimeDiff.lastTime = Date.now();
    frameTimeDiff.startBegin = Date.now();
    game();
  }
}


function platformDistance(h, heroInstance) {

  let a =
    -(2 * heroInstance.yJump) /
    (heroInstance.xJump / (2 * heroInstance.vx)) ** 2 /
    2;
  let b =
    (2 * heroInstance.yJump) / (heroInstance.xJump / (2 * heroInstance.vx));
  let delta = b ** 2 + 4 * a * h;

  return ((-b - Math.sqrt(delta)) / (2 * a)) * heroInstance.vx;
}

function level1(gridInstance, heroInstance) {
  gridInstance.grid = [];

  let d1 = platformDistance(1, heroInstance);
  let d2 = platformDistance(2, heroInstance);
  let d27 = platformDistance(2.7, heroInstance);
  let d0 = 4;
  let dMoins1 = platformDistance(-1, heroInstance);

  let pos;
  let lastPos;

  for (let k = 0; k < 70; k++) {
    gridInstance.addPlatform(k, 4);
  }

  for (let k = 70; k < 110; k++) {
    gridInstance.addPlatform(k, 4);
    if ((k - 70) % 8 === 0) {
      pos = gridInstance.addPeak(k, 5, "up");
    }
  }
  lastPos = pos;

  for (let k = 0; k < 8; k++) {
    pos = gridInstance.addPlatform(Math.floor(30 + lastPos[0]) + d0 * k, 5);
  }

  for (let k = 110; k < 188; k++) {
    pos = gridInstance.addPlatform(k, 4);
  }

  lastPos = pos;
  for (let k = 0; k < 3; k++) {
    pos = gridInstance.addPlatform(lastPos[0] + d1 * (k + 1), 5 + k);
  }
  pos = gridInstance.addPlatform(pos[0] + d0, 7);
  pos = gridInstance.addPlatform(pos[0] + dMoins1, 6);
  pos = gridInstance.addPlatform(pos[0] + d0, 6);
  pos = gridInstance.addPlatform(pos[0] + d1, 7);
  pos = gridInstance.addPlatform(pos[0] + d0, 7);
  pos = gridInstance.addPlatform(pos[0] + 2.5, 6);
  pos = gridInstance.addPlatform(pos[0] + 2.5, 5);
  for (let k = 0; k < 5; k++) {
    pos = gridInstance.addPlatform(pos[0] + d1, 6 + k);
  }
  for (let k = 0; k < 5; k++) {
    pos = gridInstance.addPlatform(pos[0] + dMoins1, 10 - k);
  }
  for (let k = 0; k < 15; k++) {
    pos = gridInstance.addPlatform(pos[0] + d1, 6 + k);
  }

  lastPos = pos;
  for (
    let k = Math.floor(lastPos[0]) + 4;
    k < Math.floor(lastPos[0]) + 24;
    k++
  ) {
    pos = gridInstance.addPlatform(k, 4);
  }

  lastPos = pos;
  for (let k = Math.floor(lastPos[0]); k < Math.floor(lastPos[0]) + 10; k++) {
    pos = gridInstance.addPlatform(k, 4);
    pos = gridInstance.addPlatform(k, 7.6);
    pos = gridInstance.addPeak(k, 6.5, "down");
  }

  pos = gridInstance.addPeak(pos[0] + 4, 5, "left");
  pos = gridInstance.addPlatform(pos[0] + 1, 5);
  pos = gridInstance.addPlatform(Math.floor(pos[0]) + 1, 4);
  pos = gridInstance.addPeak(pos[0], 5, "right");

  lastPos = pos;
  for (let k = Math.floor(lastPos[0]); k < Math.floor(lastPos[0]) + 10; k++) {
    pos = gridInstance.addPlatform(k, 4);
  }

  lastPos = pos;
  for (let k = 1; k < 5; k++) {
    pos = gridInstance.addPlatform(lastPos[0] + k * d2, 4 + 2 * k);
  }
  lastPos = pos;
  for (let k = 1; k < 5; k++) {
    pos = gridInstance.addPlatform(lastPos[0] + k * d27, lastPos[1] + 2.7 * k);
  }

  pos = gridInstance.addPlatform(pos[0] + d0, pos[1]);
  pos = gridInstance.addPlatform(pos[0] + d0, pos[1]);
  pos = gridInstance.addPlatform(pos[0] + 2.5, pos[1] - 1);
  pos = gridInstance.addPlatform(pos[0] + d1, pos[1] + 1);
  pos = gridInstance.addPlatform(pos[0] + dMoins1, pos[1] - 1);
  pos = gridInstance.addPlatform(pos[0] + 1, pos[1]);
  pos = gridInstance.addPlatform(pos[0] + 1, pos[1]);
  pos = gridInstance.addPlatform(pos[0] + 1, pos[1]);
  pos = gridInstance.addPlatform(pos[0] + 1, pos[1]);
  pos = gridInstance.addPeak(pos[0], pos[1] + 1, "up");
  pos = gridInstance.addPlatform(pos[0] + dMoins1 - 1, pos[1] - 2);
  pos = gridInstance.addPlatform(pos[0] + 2.5, pos[1] - 1);
  pos = gridInstance.addPlatform(pos[0] + 2.5, pos[1] - 1);
  pos = gridInstance.addPeak(pos[0] - 1, pos[1] + 3, "left");
  pos = gridInstance.addPlatform(pos[0] + 1 + dMoins1, pos[1] - 4);
  pos = gridInstance.addPlatform(pos[0] + dMoins1, pos[1] - 1);
  pos = gridInstance.addPlatform(pos[0] + 2.5, pos[1] - 1);
  pos = gridInstance.addPlatform(pos[0] + dMoins1, pos[1] - 1);

  for (let k = 0; k < 12; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, pos[1] - Math.abs(Math.sin(k)));
  }

  lastPos = pos;
  for (let k = lastPos[0] + 1; k < Math.floor(lastPos[0]) + 4; k++) {
    pos = gridInstance.addPlatform(k, lastPos[1]);
  }

  //gameParameters.initial[0] = [pos[0]-6, pos[1] + 10] ;
  pos = gridInstance.addPlatform(pos[0] + d0, pos[1]);

  let heightVariation;
  heightVariation = 2;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = 2.7;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = -4;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = 1;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = 1.5;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = -1;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = 2.3;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = 1;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = -4;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  heightVariation = -4;
  pos = gridInstance.addPlatform(
    pos[0] + platformDistance(heightVariation, heroInstance),
    pos[1] + heightVariation
  );
  pos = gridInstance.addPlatform(pos[0] + d0, 4);

  lastPos = pos;
  for (let k = 1; k < 20; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }

  pos = gridInstance.addPlatform(pos[0] + 1, 4);
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 5; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 5; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 6; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 4; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 4; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 8; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 1; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 7; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 1; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 3; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 1; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPeak(pos[0], 5, "up");
  for (let k = 0; k < 11; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPlatform(pos[0] + 1, 4.5);
  pos = gridInstance.addPeak(pos[0], 5.5, "up");
  for (let k = 0; k < 5; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPlatform(pos[0] + 1, 4.5);
  pos = gridInstance.addPeak(pos[0], 5.5, "up");
  for (let k = 0; k < 5; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPlatform(pos[0] + 1, 4.5);
  pos = gridInstance.addPeak(pos[0], 5.5, "up");
  for (let k = 0; k < 3; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPlatform(pos[0] + 1, 4.5);
  pos = gridInstance.addPeak(pos[0], 5.5, "up");
  for (let k = 0; k < 8; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPlatform(pos[0] + 1, 5);
  pos = gridInstance.addPeak(pos[0], 6, "up");

  for (let k = 0; k < 8; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }
  pos = gridInstance.addPlatform(pos[0] + d0, 4);
  pos = gridInstance.addPlatform(pos[0] + d0, 4);
  pos = gridInstance.addPlatform(pos[0] + 1, 4);
  pos = gridInstance.addPlatform(pos[0] + 1, 4);
  pos = gridInstance.addPlatform(pos[0] + 1, 4);
  pos = gridInstance.addPlatform(pos[0] + dMoins1, 3);
  pos = gridInstance.addPlatform(pos[0] + d0, 3);
  pos = gridInstance.addPlatform(pos[0] + dMoins1, 2);
  pos = gridInstance.addPlatform(pos[0] + 1, 2);
  pos = gridInstance.addPlatform(pos[0] + 1, 2);
  pos = gridInstance.addPlatform(pos[0] + 1, 2);
  lastPos = pos;
  for (let k = 0; k < 16; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 2);
  }
  pos = lastPos;
  for (let k = 0; k < 16; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4.3);
  }

  for (let l = 2; l < 16; l += 2) {
    lastPos = pos;
    for (let k = 0; k < 3; k++) {
      pos = gridInstance.addPlatform(pos[0] + 1, l);
    }
    pos = gridInstance.addPlatform(pos[0] + 1, l + 2);
    pos = gridInstance.addPlatform(pos[0], l);
    pos = gridInstance.addPlatform(pos[0], l + 1);
    pos = gridInstance.addPlatform(pos[0], l + 2);
    pos = lastPos;
    if (l == 2) {
      pos = gridInstance.addPlatform(pos[0], 5.3 + (l - 2));
    }
    pos = gridInstance.addPlatform(pos[0], 5.3 + (l - 2) + 1);
    for (let k = 0; k < 3; k++) {
      pos = gridInstance.addPlatform(pos[0], 5.3 + (l - 2) + k + 1);
    }
    for (let k = 0; k < 3; k++) {
      pos = gridInstance.addPlatform(pos[0] + 1, 6.3 + l);
      pos = gridInstance.addPeak(pos[0], 5.3 + l, "down");
    }
  }

  for (let k = 1; k < 20; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 16);
    pos = gridInstance.addPlatform(pos[0], 20.3);
    pos = gridInstance.addPeak(pos[0], 19.3, "down");
  }

  for (let k = 1; k < 40; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 16);
    if (k % 5 === 0) {
      pos = gridInstance.addPeak(pos[0], 17.5, "left");
    }
  }
  for (let k = 1; k < 40; k++) {
    if (k % 5 === 0) {
      pos = gridInstance.addPlatform(pos[0] + d0, 16);
    }
  }
  for (let k = 1; k < 15; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 16);
  }
  for (let k = 1; k < 40; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 16);
    if (k % 5 === 0) {
      lastPos = pos;
      pos = gridInstance.addPeak(pos[0] + Math.sin(k), 17.5, "left");
      pos = gridInstance.addPeak(pos[0] + 0.5, 17.5, "left");
      pos = lastPos;
    }
  }
  for (let k = 1; k < 10; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 16);
  }
  for (let k = 0; k < 20; k++) {
    pos = gridInstance.addPlatform(pos[0] + d1, 16 + k);
  }
  lastPos = pos;
  pos = gridInstance.addPlatform(pos[0] + d1 + 1, pos[1] + 2);
  pos = gridInstance.addPeak(pos[0] - 1, pos[1], "left");
  pos = lastPos;

  for (let k = 0; k < 8; k++) {
    pos = gridInstance.addPlatform(pos[0] + 2.5, pos[1] - 1);
  }
  for (let k = 0; k < 60; k++) {
    pos = gridInstance.addPlatform(pos[0] + d1, pos[1] + 1);
  }

  pos = gridInstance.addPlatform(pos[0] + 8, pos[1] - 22);

  for (let k = 0; k < 7; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, pos[1]);
  }

  pos = gridInstance.addPlatform(pos[0] + 8, pos[1] - 22);

  for (let k = 0; k < 7; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, pos[1]);
  }
  pos = gridInstance.addPeak(pos[0], pos[1] + 1, "up");

  for (let k = 1; k < 15; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }

  for (let k = 0; k < 40; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
    if (k % 4 === 0) {
      pos = gridInstance.addPeak(pos[0], 5, "up");
    }
  }

  for (let k = 0; k < 40; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
    if (k % 8 === 0) {
      pos = gridInstance.addPeak(pos[0], 5, "up");
    }
  }

  for (let k = 0; k < 80; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
    if (k % 16 === 0) {
      pos = gridInstance.addPeak(pos[0], 5, "up");
    }
  }

  lastPos = pos;
  for (let k = 1; k < 80; k++) {
    pos = gridInstance.addPlatform(pos[0] + 1, 4);
  }

  let endingInstance = new ending(Math.floor(lastPos[0]) + 33);
  gridInstance.addEnding(endingInstance);
}


const keys = {};

const frameTimeDiff = {
  lastTime: 0,

  dt: 0,

  startBegin: 0,

  endingBegin: 0,

  lastCheckPoint: 0,

};


const checkPointValue = {};


const backGroundPositionSauv = { save: 0 };

let checkPointCounter = 0;


const gameParameters = {
  initial: [[6, 5.5], [1, 0], 10, 0, 4, 3, 0, 0, false],
};


let heroInstance;


let drawingInstance;


let gridInstance;

function AtLoad() {

  heroInstance = new hero(
    gameParameters.initial[0],
    gameParameters.initial[1],
    gameParameters.initial[2],
    gameParameters.initial[3],
    gameParameters.initial[4],
    gameParameters.initial[5],
    gameParameters.initial[6],
    gameParameters.initial[7],
    gameParameters.initial[8]
  );
  // construct hero
  drawingInstance = new drawing(heroInstance);
  gridInstance = new grid();
  level1(gridInstance, heroInstance);

  checkPointValue.checkpoint = new checkPoint(heroInstance);

  drawingInstance.drawMovingtext(
    "Press Space to begin",
    0,
    0,
    0,
    0,
    80,
    [0, 16],
    [50, 16],
    heroInstance
  );
  setTimeout(function () {
    restart(gameParameters.initial);
  }, 1000);
}

window.onload = AtLoad;

window.addEventListener("keydown", keyEventHandler);
window.addEventListener("keyup", keyEventHandler);
