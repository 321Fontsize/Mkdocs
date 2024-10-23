document.addEventListener('DOMContentLoaded', function () {
    // 获取footer的meta区域
    var footerMetaInner = document.querySelector('.md-footer-meta__inner.md-grid');
  
    // 创建新的span元素
    var busuanziContainer = document.createElement('span');
    busuanziContainer.id = 'busuanzi_container_page_pv';
    busuanziContainer.innerHTML = '总阅读量<span id="busuanzi_value_page_pv"></span>次';
  
    // 将新的span元素添加到footerMetaInner中
    if (footerMetaInner) {
      footerMetaInner.appendChild(busuanziContainer);
    }
  });