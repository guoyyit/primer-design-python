# 本地 DNA 引物设计与特异性评估服务

给定模板序列与设计参数，服务在本地完成候选引物搜索、理化性质计算、非目标结合风险
评估，并返回**可复现**的排序结果。纯 Python 实现（依赖 FastAPI），无需外部二进制
或联网数据库。

## 快速开始

```bash
uv sync --dev
uv run uvicorn main:app                # 启动服务，默认数据目录 ./.primer_data
PRIMER_DATA_DIR=/var/lib/primer uv run uvicorn main:app --port 8000
uv run pytest                           # 运行测试
```

Python 库直接使用：

```python
from primer_design.service import PrimerDesignService
from primer_design.models import DesignParams

svc = PrimerDesignService("./data")     # 传 None 可完全禁用持久化
svc.start()
task = svc.submit_design(DesignParams(template="ACGT...", circular=False))
print(svc.get_task(task.task_id).status, svc.get_task(task.task_id).progress)
```

## 设计流程

1. **候选生成**：在目标区间 `[target_start, target_end)` 内枚举所有长度位于
   `[min_primer_len, max_primer_len]` 的正向窗口；反向引物为窗口的反向互补。
   线性模板坐标超出末端的窗口被截断/丢弃；**环状模板**允许窗口跨 0 点环绕，
   扩增子长度按环绕距离 `(rev_end - fwd_start) mod n` 计算。
2. **单条过滤**：GC 含量、退火温度、自身二聚体、发夹必须全部落在阈值内。
3. **配对过滤**：正反向引物方向正确（反向引物位于正向引物下游），产物长度位于
   `[min_product_len, max_product_len]`，异二聚体稳定性不超过 `max_dimer_dg`。
4. **结构排序**：对通过过滤的配对计算综合惩罚分（见下），先取结构短名单
   （`max_candidates × 3`）。
5. **特异性评估**：对短名单逐个在参考库当前快照上扫描非目标结合，把命中风险
   计入最终惩罚后重新排序，返回前 `max_candidates` 个。

## 理化性质与判定规则

| 指标 | 方法 | 判定/阈值字段 |
| --- | --- | --- |
| 退火温度 Tm | SantaLucia (1998) 统一最近邻参数，含起始项与末端 AT 罚分；无盐校正，引物浓度 `tm_concentration_nm` 按非自互补 `Ct/4` 近似 | `min_tm` / `max_tm` |
| GC 含量 | `(G+C)/L` | `min_gc` / `max_gc` |
| 自身二聚体 | 枚举所有反平行末端重叠对齐，最近邻堆积能求和，取最负 ΔG（kcal/mol，37 ℃） | `max_dimer_dg` |
| 异二聚体 | 同一对齐枚举应用于正反向两条引物 | `max_dimer_dg` |
| 发夹 | 枚举所有外部茎配对，要求茎长 ≥ 3 bp、环 ≥ 3 nt，取最负茎 ΔG + 0.1×环长 | `max_hairpin_dg` |
| 连续重复 | 最长单碱基同聚物；最长二核苷酸串联重复（如 ATATAT） | 排序惩罚（默认 > 4 / > 3） |
| 末端回文 | 两端各 4 碱基互为反向互补 | 排序惩罚（固定罚分） |
| 3′ GC 夹 | 末 5 个碱基中 G/C 数 | 排序惩罚（> 3 时） |

ΔG 约定：**越负越稳定，越不希望出现**。过滤语义为
`ΔG >= max_dimer_dg` 才保留（例如阈值 `-6.0` 表示比 −6.0 kcal/mol 更稳定的
结构被拒绝；阈值设为正数会拒绝一切配对）。

### 排序权重（`primer_design/scoring.py`）

总分 = 以下各项之和，**分数越低越好**。同一总分按
`(分数, 产物长度, 正向序列, 反向序列)` 决胜，保证相同输入（含 `seed`）逐字节一致：

| 惩罚项 | 权重 | 说明 |
| --- | --- | --- |
| `tm_distance` | 1.0 / ℃ | 每条引物相对 Tm 窗口中心的距离 |
| `gc_distance` | 8.0 / 单位 | GC 含量偏离 0.50 |
| `self_dimer` | 1.0 / kcal/mol | 低于 −3.0 地板的部分 |
| `hairpin` | 1.0 / kcal/mol | 低于 −1.0 地板的部分 |
| `homopolymer` | 1.5 / 碱基 | 超过 4 的部分 |
| `dinucleotide_run` | 1.0 / 重复 | 超过 3 个重复的部分 |
| `terminal_palindrome` | 4.0 固定 | 命中即罚 |
| `gc_clamp` | 1.0 / 碱基 | 末 5 碱基中超过 3 个 G/C 的部分 |
| `tm_delta` | 1.5 / ℃ | 正反向 Tm 差 |
| `product_len` | 0.02 / 碱基 | 产物长度偏离请求区间中心 |
| `heterodimer` | 1.5 / kcal/mol | 低于 −3.0 地板的部分 |
| 非目标命中 | 每条 ≥ 3.0 | 基础罚分 + 错配越少越重 + 3′ 连续匹配越长越重 |

每条返回结果都带 `penalty_detail`，可直接核对各项贡献。

## 非目标结合（特异性）

* 对每条引物，在参考序列（正向与反向互补两条链）上做**逐窗口全比对**：
  窗口错配数 ≤ `max_mismatches`，且引物 **3′ 端连续匹配数** ≥
  `min_three_prime_match` 才记为一次命中。
* 长参考序列按 `chunk_size` 分块扫描（相邻块重叠 `L−1` 碱基，跨块窗口不丢），
  额外内存恒定，进度按已扫描碱量上报。
* 命中结果按 `(ref_id, strand, position)` 排序；每条引物每个参考最多保留 100 条，
  避免极端重复库撑爆结果。
* **快照一致性**：参考库以不可变元组 + 版本号保存。更新（新增/替换/清空）先构造
  全新快照再原子替换；一次评估开始时取一个快照，期间发生的库更新对其不可见，
  两条引物看到的也一定是同一版本，不会读到新旧混合状态。

## 任务、取消与持久化

* 任务经后台 worker 串行（默认 `max_workers=1`）执行，进度区间：
  `0–0.2` 候选生成、`0.2–0.6` 配对枚举、`0.6–1.0` 特异性评估。
* `cancel` 对排队任务立即标记取消；运行中任务在下一个取消检查点（枚举循环与
  每个参考分块）抛出 `CancelledError`，结果与中间产物被丢弃，不产生部分写入。
* 数据目录内所有写盘均为 **临时文件 + fsync + `os.replace`** 原子操作。
  启动时清理 `*.partial` 残留；取消或崩溃不会留下临时索引。
* 重启恢复：已完成任务（含完整结果）直接可查询；排队/运行中的任务被重置为
  queued、进度归零并重新执行；**中间结果从不落盘**，因此不可能被当作最终结果返回。
* 参考库与任务分文件持久化（`library.json` / `tasks.json`），损坏的单条任务
  记录会被跳过而不影响其余记录恢复。

## HTTP API

| 方法/路径 | 说明 |
| --- | --- |
| `POST /designs` | 提交设计任务（202），立即校验序列与参数区间 |
| `GET /designs` | 列出所有任务 |
| `GET /designs/{id}` | 查询状态、进度与完成后的排序结果 |
| `POST /designs/{id}/cancel` | 取消任务（终态返回 409） |
| `GET /references` | 库版本与条数 |
| `POST /references` | 增量新增（重复 id 返回 409） |
| `PUT /references` | 全量替换 |
| `DELETE /references` | 清空 |
| `GET /healthz` | 健康检查 |

## 可复现性

排序是纯函数：无随机抽样，候选按 `(position, length)` 枚举，命中按坐标排序，
最终按上述全序键决胜。相同模板、参数与相同参考库快照下，输出逐条一致；
`seed` 参数被保留在参数模型中（当前确定性算法不直接使用），供未来引入随机搜索
时保持承诺。

## 限制（刻意的工程取舍）

* 结构能量是**近似最近邻堆积评分**，不是 mfold/UNAFold 完整自由能：不建模
  凸环、多分支环、盐浓度对 Tm/ΔG 的影响，也不使用实验外推的发夹参数。
* 特异性是逐窗口精确比对，时间复杂度 O(L·R)：对百万级参考库属于"能跑完但不
  实时"。实现刻意选择正确性与内存可控，而非构建大型 k-mer 索引。
* 仅接受 A/C/G/T 作为引物窗口；模板支持 IUPAC 歧义碱基，但含歧义碱基的窗口
  不会成为候选。
* 配对枚举有 200000 对上限以防极端参数失控；命中每参考每引物最多 100 条。
* worker 默认单线程；大批量请求建议分块提交并通过进度接口轮询。
