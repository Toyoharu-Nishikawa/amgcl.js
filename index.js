import AmgCL from "./amgcl.js"

export const version = ()=>"AMGCL version 1.4.3"

export const AMGCL = async () => {
  const amgcl = await AmgCL()

  const func = {
    init(){
      const amg = new AMG(amgcl)
      return amg 
    },
    fromMatrix(matrix){
      const amg = new AMG(amgcl)
      amg.setMatrix(matrix)
      return amg 
    },
    fromCRS(m,n,nnz,a,asub,xa){
      const amg = new AMG(amgcl)
      amg.setCRS(m,n,nnz,a,asub,xa)
      return amg
    }
  }
  return func 
}

const AMG = class{
  constructor(amgcl){
    this.amgcl = amgcl
    this.rows = null
    this.nnz = null
    this.ptr = null 
    this.col = null
    this.val = null
    this.rhs = null
    this.iniX= null	  
    this.maxIteration = null
    this.tolerance = null

    this.matrix = null 
  }
  get CRS(){
    const rows = this.rows
    const nnz  = this.nnz 
    const ptr  = this.ptr 
    const col  = this.col 
    const val  = this.val 
 
    const obj  = {rows, nnz, ptr, col, val}
    return obj
  }
  set CRS(obj){
    this.rows = obj.rows
    this.nnz  = obj.nnz 
    this.ptr  = obj.ptr 
    this.col  = obj.col 
    this.val  = obj.val 
  }
  clear(){
    this.rows = null
    this.nnz = null
    this.ptr = null
    this.col = null
    this.val = null
    this.rhs = null
    this.iniX=null  
    this.maxIteration = null
    this.tolerance = null
    this.matrix = null 

    return this
  }
  release(){
    this.superlu = null
    this.rows = null
    this.nnz = null
    this.ptr = null
    this.col = null
    this.val = null
    this.rhs = null
    this.iniX=null  
    this.maxIteration = null
    this.tolerance = null
    this.matrix = null 

    return this
  }
  setMatrix(A){
    this.matrix = A
    this.setCRSFromMatrix(A)

    return this
  }
  setCRS(rows,nnz,ptr,col,val){
    this.rows = obj.rows
    this.nnz  = obj.nnz 
    this.ptr  = obj.ptr 
    this.col  = obj.col 
    this.val  = obj.val 


    return this
  }
  setRHS(B){
    if(!Array.isArray(B)){
      return
    }
    if(Array.isArray(B[0])){
      this.nrhs = B.length
      this.rhs = this.convert2Dto1D(B)
    }
    else{
      this.rhs = B
      this.nrhs = 1
    }

    return this
  }
  solve(rhs,tolerance,maxIteration,initialX){
    const amgcl  = this.amgcl
    const rows = this.rows
    const nnz  = this.nnz 
    const ptr  = this.ptr 
    const col  = this.col 
    const val  = this.val 
 
    this.rhs = rhs
    const iniX = Array.isArray(initialX)? initialX :[...Array(rows)].fill(0)
    this.iniX = iniX
    this.tolerance = tolerance
    this.maxIteration = maxIteration

    const res = linSolve(amgcl, rows, nnz,ptr,col,val,rhs,tolerance,maxIteration,iniX)

    return res 
  }
  setCRSFromMatrix(A){ // Compressed Row Storage
    const rows = A.length    
    
    const val = []
    const col = []
    const ptr = [0]
    let count =0
    for(let i=0;i<rows;i++){
      for(let j=0;j<rows;j++){
        const aij = A[i][j]
        if(aij===0){continue}
        val.push(aij)
        col.push(j)
        ++count
      }
      ptr.push(count)
    }
    const nnz = val.length
    this.rows = rows
    this.nnz = nnz
    this.val = val
    this.col = col
    this.ptr = ptr

    return this
  }
}

const linSolve = (amgcl,rows, nnz,ptr_,col_,val_,rhs_,tolerance,maxIteration,iniX_) => {
  const F64Byte =  8 
  const I32Byte =  4 

  const ptr   = new Int32Array([...ptr_])
  const col   = new Int32Array([...col_])
  const val   = new Float64Array([...val_])
  const rhs   = new Float64Array([...rhs_])
  const iniX  = new Float64Array([...iniX_])
  const X     = new Float64Array(rows)
  const err   = new Float64Array(1)
  const itr   = new Int32Array(1)

  const ptrP    = amgcl._malloc((rows+1) * I32Byte)
  const colP    = amgcl._malloc(nnz      * I32Byte)
  const valP    = amgcl._malloc(nnz      * F64Byte)
  const rhsP    = amgcl._malloc(rows     * F64Byte)
  const iniXP   = amgcl._malloc(rows     * F64Byte)
  const XP      = amgcl._malloc(rows     * F64Byte)
  const errP    = amgcl._malloc(1        * F64Byte)
  const itrP    = amgcl._malloc(1        * I32Byte)


  amgcl.HEAP32.set(ptr   , ptrP/I32Byte)
  amgcl.HEAP32.set(col   , colP/I32Byte)
  amgcl.HEAPF64.set(val  , valP/F64Byte)
  amgcl.HEAPF64.set(rhs  , rhsP/F64Byte)
  amgcl.HEAPF64.set(iniX , iniXP/F64Byte)
  amgcl.HEAPF64.set(X    , XP/F64Byte)
  amgcl.HEAPF64.set(err  , errP/F64Byte)
  amgcl.HEAP32.set(itr   , itrP/I32Byte)

  amgcl._solve(rows,nnz, ptrP, colP,valP, rhsP,iniXP,tolerance,maxIteration, XP,errP, itrP )


  const itrArray = new Int32Array(amgcl.HEAP32.buffer, itrP,1)
  const errArray = new Float64Array(amgcl.HEAPF64.buffer, errP,1)
  const xArray =   new Float64Array(amgcl.HEAPF64.buffer, XP,rows)

  const iteration = itrArray[0]
  const residualError = errArray[0]
  const x = [...xArray]
  const res = {iteration, residualError, x}

  amgcl._free(ptrP )
  amgcl._free(colP )
  amgcl._free(valP )
  amgcl._free(rhsP )
  amgcl._free(iniXP)
  amgcl._free(XP   ) 
  amgcl._free(errP ) 
  amgcl._free(itrP ) 


  return res 
}


