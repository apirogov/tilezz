#[cfg(test)]
#[allow(clippy::module_inception)]
mod unit_tests {
    use crate::cyclotomic::Units;
    use crate::cyclotomic::*;

    macro_rules! test_units {
        ($t:ty) => {{
            type Z = $t;
            let turn = Z::turn();
            // periodicity
            for k in 0..turn {
                assert_eq!(<Z as Units>::unit(k), <Z as Units>::unit((k + turn) as i8));
                assert_eq!(<Z as Units>::unit(k), <Z as Units>::unit((k - turn) as i8));
            }

            // group property: unit(a)*unit(b)=unit(a+b)
            for a in 0..turn {
                for b in 0..turn {
                    let lhs = <Z as Units>::unit(a) * <Z as Units>::unit(b);
                    let rhs = <Z as Units>::unit(((a as i16 + b as i16) % (turn as i16)) as i8);
                    assert_eq!(lhs, rhs);
                }
            }

            // ccw is unit(1)
            assert_eq!(Z::ccw(), <Z as Units>::unit(1));
            assert_eq!(Z::one(), <Z as Units>::unit(0));
        }};
    }

    #[test]
    fn test_unit_tables() {
        test_units!(ZZ4);
        test_units!(ZZ6);
        test_units!(ZZ8);
        test_units!(ZZ10);
        test_units!(ZZ12);
        test_units!(ZZ16);
        test_units!(ZZ20);
        test_units!(ZZ24);
        test_units!(ZZ30);
        test_units!(ZZ32);
        test_units!(ZZ60);
    }
}
