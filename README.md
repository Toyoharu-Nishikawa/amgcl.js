# amgcl.js


```samlple.js
import {AMGCL} from "./packages/index.js"

export let AMG = null
const main = async () => {
  const s = 19
  const u = 21
  const p = 16
  const e =  5
  const r = 18
  const l = 12

  const A = [
    [s, 0, u, u, 0],
    [l, u, 0, 0, 0],
    [0, l, p, 0, 0],
    [0, 0, 0, e, u],
    [l, l, 0, 0, r],
  ]

  AMG = await AMGCL()
  console.log("AMG", AMG)

  /** 
   * set system matrix from double array format
   */
  const amg = AMG.fromMatrix(A)
  
  /** set system matrix from CRS(compressed row strage) format
    const col = [0, 2, 3, 0, 1, 1, 2, 3, 4, 0, 1, 4]
    const nnz = 12
    const ptr = [0, 3, 5, 7, 9, 12]
    const rows = 5
    const val = [19, 21, 21, 12, 21, 12, 16, 5, 21, 12, 12, 18]

    const amg = AMG.fromCRS(rows,nnz,ptr,col,val)
   */

  console.log("amg", amg)
  console.log("CRS",amg.CRS)


  /**
   * right hand side vector
   */
  const V = [1,1,1,1,1]


  const  ans = [-0.031250,  0.065476,  0.013393,  0.062500,  0.032738]
  console.log("ans", ans)

  const tolerance = 1E-8
  const maxIteration = 1E3
  const iniX = [0,0,0,0,0]
  const res = amg.solve(V,tolerance,maxIteration,iniX) //  solve for x: Ax = V^T (B is trnasposed)

  console.log("res",res)
  const x = res.x
}
main()
```
