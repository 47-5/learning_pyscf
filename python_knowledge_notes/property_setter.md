# Python 知识点：`@property` 和 `@xxx.setter`

## 背景代码

```python
@nelec.setter
def nelec(self, neleca_nelecb):
    neleca, nelecb = neleca_nelecb
    self._nelectron = neleca + nelecb
    self.spin = neleca - nelecb
```

这里的 `@nelec.setter` 是给前面已经定义好的 `nelec` 属性添加“赋值时的行为”。

## `@property` 做什么

`@property` 可以把一个方法包装成“像属性一样访问”的东西。

例如：

```python
class Molecule:
    @property
    def nelec(self):
        return self._nelectron
```

这样外部就可以写：

```python
mol.nelec
```

而不是：

```python
mol.nelec()
```

也就是说，`@property` 主要控制“读取属性时发生什么”。

## `@nelec.setter` 做什么

如果只写了 `@property`，这个属性默认通常只能读，不能直接赋值。

比如：

```python
mol.nelec = (5, 4)
```

如果想让这个赋值动作触发自定义逻辑，就需要写 setter：

```python
class Molecule:
    @property
    def nelec(self):
        return self._nelectron

    @nelec.setter
    def nelec(self, neleca_nelecb):
        neleca, nelecb = neleca_nelecb
        self._nelectron = neleca + nelecb
        self.spin = neleca - nelecb
```

现在执行：

```python
mol.nelec = (5, 4)
```

实际会调用：

```python
Molecule.nelec.fset(mol, (5, 4))
```

从效果上看，相当于：

```python
neleca, nelecb = (5, 4)
mol._nelectron = 9
mol.spin = 1
```

## 为什么这个例子适合用 setter

`nelec` 表面上像一个普通属性，但它背后代表的是两个相关量：

- 总电子数：`_nelectron = neleca + nelecb`
- 自旋：`spin = neleca - nelecb`

所以写：

```python
mol.nelec = (neleca, nelecb)
```

不仅仅是保存一个二元组，而是顺便更新对象内部的两个状态。

这就是 setter 常见的用途：让“属性赋值”不只是简单赋值，而是带有校验、转换、同步其他字段等逻辑。

## 一句话记忆

`@property` 管读取，`@xxx.setter` 管赋值。

对应关系是：

```python
obj.xxx        # 调用 @property 那个函数
obj.xxx = val  # 调用 @xxx.setter 那个函数
```

## 一个最小例子

```python
class Circle:
    def __init__(self, radius):
        self._radius = radius

    @property
    def radius(self):
        return self._radius

    @radius.setter
    def radius(self, value):
        if value <= 0:
            raise ValueError("radius must be positive")
        self._radius = value
```

使用：

```python
c = Circle(2)
print(c.radius)  # 读取，调用 @property

c.radius = 3     # 赋值，调用 @radius.setter
c.radius = -1    # 报错，因为 setter 里做了检查
```

## 小心点

setter 函数的名字必须和 property 的名字一样。

正确：

```python
@property
def nelec(self):
    ...

@nelec.setter
def nelec(self, value):
    ...
```

不要写成：

```python
@nelec.setter
def set_nelec(self, value):
    ...
```

否则类里的属性名字关系会变得混乱。
